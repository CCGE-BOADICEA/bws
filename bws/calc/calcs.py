"""
Risk and mutation probability calculations, for details
see https://github.com/CCGE-BOADICEA/boadicea/wiki/Cancer-Risk-Calculations

© 2022 Cambridge University
SPDX-FileCopyrightText: 2022 Cambridge University
SPDX-License-Identifier: GPL-3.0-or-later
"""
from bws import pedigree
from bws.exceptions import TimeOutException, ModelError
from collections import OrderedDict
from django.conf import settings
from django.http.request import HttpRequest
from rest_framework.exceptions import ValidationError
from rest_framework.request import Request
from subprocess import Popen, PIPE, TimeoutExpired
import logging
import os
import re
import resource
import tempfile
import time
from bws.calc.model import ModelParams, ModelOpts
from bws.calc.risks import Risk, RemainingLifetimeBaselineRisk, RiskBaseline
from bws.pedigree import Pedigree


logger = logging.getLogger(__name__)
REGEX_ALPHANUM_COMMAS = re.compile("^([\\w,]+)$")


class Predictions(object):

    def __init__(self, pedi, model_params=ModelParams(),
                 risk_factor_code=0, hgt=-1, mdensity=None, prs=None, cwd=None, request=Request(HttpRequest()),
                 run_risks=True, model_settings=settings.BC_MODEL, calcs=None):
        """
        Run cancer risk and mutation probability prediction calculations.
        @param pedi: L{Pedigree} used in prediction calculations
        @keyword model_params: L{ModelParams}, model parameters
        @keyword risk_factor_code: risk factor code
        @keyword hgt: height
        @keyword mdensity: mammographic density  
        @keyword prs: polygenic risk alpha & beta values calculated from VCF file
        @keyword cwd: working directory
        @keyword request: HTTP request
        @keyword run_risks: run risk calculations, default True
        @keyword model_settings: cancer model settings
        @keyword calcs: list of calculations to run, e.g. ['carrier_probs', 'remaining_lifetime']
        """
        assert isinstance(pedi, Pedigree), "%r is not a Pedigree" % pedi
        assert isinstance(model_params, ModelParams), "%r is not a ModelParams" % model_params
        self.pedi = pedi
        self.model_params = model_params
        self.request = request
        self.cwd = cwd
        self.risk_factor_code = risk_factor_code
        self.hgt = hgt
        self.mdensity = mdensity
        self.prs = prs
        self.model_settings = model_settings
        self.calcs = self.model_settings['CALCS'] if calcs is None else calcs

        for c in self.calcs:
            if c not in settings.ALLOWED_CALCS:     # check calculations are in the allowed list
                raise ValidationError("Unknown calculation requested: "+c)

        if cwd is None:
            self.cwd = tempfile.mkdtemp(prefix=str(request.user)+"_", dir="/tmp")
        if isinstance(risk_factor_code, int):
            self.risk_factor_code = str(risk_factor_code)
        if run_risks:
            self._run_risks()

    def is_calculate(self, calc):
        '''
        Determine if a calculation is to be run.
        @param calc: calculation name, e.g. 'carrier_probs', 'remaining_lifetime'
        '''
        return True if len(self.calcs) == 0 else (calc in self.calcs)

    def _run_risk(self, risk, model_opts):
        """
        Calculate the risk and return the parsed output as a list.
        @return: list of risks for each age
        """
        p = risk.get_pedigree()
        pf = p.write_pedigree_file(risk_factor_code=risk.get_risk_factor_code(),
                                   hgt=risk.get_hgt(),
                                   mdensity=risk.get_md(),
                                   prs=risk.get_prs(),
                                   filepath=os.path.join(self.cwd, risk.type()+"_risk.ped"),
                                   model_settings=self.model_settings)
        bf = p.write_batch_file(pf,
                                filepath=os.path.join(self.cwd, risk.type()+"_risk.bat"),
                                model_settings=self.model_settings,
                                calc_ages=risk.risk_age)
        paramf = p.write_param_file(filepath=os.path.join(self.cwd, risk.type()+"_risk.params"),
                                    model_settings=self.model_settings,
                                    mutation_freq=risk.get_mutation_frequency(),
                                    isashk=self.model_params.isashk,
                                    sensitivity=self.model_params.mutation_sensitivity)
        risks = Predictions.run(self.request, bf,
                                model_opts=model_opts,
                                model_params=self.model_params,
                                param_file=paramf,
                                cwd=self.cwd,
                                niceness=self.niceness, name=risk.get_name(),
                                model=self.model_settings)
        return self._parse_risks_output(risks, model_opts)

    def _run_risks(self):
        ''' Run risk and mutation probability calculations '''
        self.version = Predictions._get_version(model=self.model_settings, cwd=self.cwd)
        self.niceness = Predictions._get_niceness(self.pedi)
        start = time.time()
        model_opts = ModelOpts.factory(self)

        rl, rr, ry, rj, mp = self._run_risk(Risk(self), model_opts)
        if rl is not None:
            self.lifetime_cancer_risk = rl
        if rr is not None:
            self.cancer_risks = rr
        if ry is not None:
            self.ten_yr_cancer_risk = ry
        if hasattr(self, "lifetime_cancer_risk") and rj is not None:
            if len(rj) > 0:
                self.ten_yr_nhs_protocol = rj
            elif ry is not None:
                self.ten_yr_nhs_protocol = ry
        if mp is not None:
            self.mutation_probabilties = mp

        # remaining lifetime baseline
        if rr is not None:
            _rl, rr, _ry, _rj, _mp = self._run_risk(RemainingLifetimeBaselineRisk(self), ModelOpts(probs=False,
                                                                                            rj=False, rl=False,
                                                                                            rr=True, ry=False))
            if rr is not None:
                self.baseline_cancer_risks = rr

            # baseline lifetime cancer risk and baseline 10-year cancer risk
            if hasattr(self, "lifetime_cancer_risk") or hasattr(self, "ten_yr_cancer_risk"):
                rl, _rr, ry, _rj, _mp = self._run_risk(RiskBaseline(self), ModelOpts(probs=False, rj=False, rr=False,
                                                                              rl=hasattr(self, "lifetime_cancer_risk"),
                                                                              ry=hasattr(self, "ten_yr_cancer_risk")))
                if rl is not None:
                    self.baseline_lifetime_cancer_risk = rl
                if ry is not None:
                    self.baseline_ten_yr_cancer_risk = ry

        name = str(self.model_settings.get('NAME', ""))
        logger.info(
            f"{name} CALCULATIONS: user={self.request.user.id}; "
            f"elapsed time={time.time() - start}; "
            f"pedigree size={len(self.pedi.people)}; "
            f"version={getattr(self, 'version', 'N/A')}")

    def _parse_risks_output(self, risks, model_opts):
        """
        Parse computed cancer risk results.
        @param risks: cancer risks text from fortran output
        @param model_opts: cancer model options
        @return: rl, lifetime cancer risk calculation (proband's age set at 20y, censor age set at 80y)
                 rr, cancer risk calculation is performed (proband's age and censor age set in input files)
                 ry, 10-yr cancer risk calculation (proband's age set at 40y, censor age set at 50y)
                 rj, 10-yr cancer risk calculation (NHS protocol for young women at high risk)
                 mp, mutation carrier probabilities
        """
        lines = risks.split(sep="\n")
        rr_arr = [] if model_opts.rr else None      # remaining lifetime risk
        rl_arr = [] if model_opts.rl else None      # lifetime risk
        ry_arr = [] if model_opts.ry else None      # 10 yr risk (40-50y)
        rj_arr = [] if model_opts.rj else None      # 10 yr risk (NHS protocol)

        rr, rl, ry, rj, mp = False, False, False, False, False
        mp_lines = ""
        model_settings = self.model_settings
        ctype = "breast" if model_settings['NAME'] == 'BC' else "ovarian"

        for _idx, line in enumerate(lines):
            if line.startswith('##'):
                rr, rl, ry, rj, mp = False, False, False, False, False
            elif 'Age' in line or pedigree.BLANK_LINE.match(line):
                continue

            if rr or rl or ry or rj:
                prts = line.split(sep=",")
                plen = len(prts)
                v = float(prts[plen-1])
                od = OrderedDict([
                        ("age", int(prts[plen-2])),
                        (ctype+" cancer risk", {"decimal": v, "percent": round(v*100, 1)})
                    ])

                if rr:
                    rr_arr.append(od)
                elif rl:
                    rl_arr.append(od)
                elif ry and ry_arr is not None:
                    ry_arr.append(od)
                elif rj:
                    rj_arr.append(od)
            elif mp:
                mp_lines += line+"\n"

            if line.startswith('## REMAINING LIFETIME RISK'):
                rr = True
            elif line.startswith('## LIFETIME RISK'):
                rl = True
            elif line.startswith('## 10-YEAR RISK') and "NHS PROTOCOL" not in line:
                ry = True
            elif line.startswith('## 10-YEAR RISK') and "NHS PROTOCOL" in line:
                rj = True
            elif line.startswith('## PROBABILITIES'):
                mp = True

        mp_arr = Predictions._parse_probs_output(mp_lines, model_settings) if model_opts.probs else None
        return rl_arr, rr_arr, ry_arr, rj_arr, mp_arr

    @classmethod
    def _parse_probs_output(cls, probs, model_settings):
        """
        Parse computed mutation carrier probability results.
        @param probs: mutation probability text from fortran output
        @param model_settings: cancer model settings
        @return: list of containing dictionaries of the mutaion probability results
        """
        probs_arr = []
        gene_columns = []

        for _idx, line in enumerate(probs.splitlines()):
            if REGEX_ALPHANUM_COMMAS.match(line):
                gene_columns = line.strip().split(sep=",")
            elif not line.startswith('#'):
                parts = line.strip().split(sep=",")

                probs_arr.append({"no mutation": {"decimal": float(parts[0]),
                                                  "percent": round(float(parts[0])*100, 2)}})
                for i, gene in enumerate(model_settings['GENES'], 1):   # 1-based loop
                    assert gene == gene_columns[i], "MUTATION CARRIER PROBABILITY - RESULTS COLUMN MISMATCH FOUND"
                    probs_arr.append({gene:
                                      {"decimal": float(parts[i]),
                                       "percent": round(float(parts[i])*100, 2)}})
        return probs_arr

    @classmethod
    def _get_niceness(cls, pedi, factor=15):
        """
        Based on the pedigree size and structure get priority setting. A positive nice value means
        lower priority. '0' means priority is not adjusted. The niceness scale is from -20 to 19.
        @param pedi: L{Pedigree} to return a priority for.
        @param factor: used to calculate a niceness value
        @return: niceness: priority
        """
        (siblings, _siblings_same_yob) = pedi.get_siblings(pedi.get_target())
        nsiblings = len(siblings)
        if nsiblings > 0:
            return nsiblings if nsiblings < 19 else 19

        pedigree_size = len(pedi.people)
        niceness = int(pedigree_size/factor)
        if niceness > 19:
            niceness = 19
        return niceness

    @classmethod
    def _get_version(cls, model=settings.BC_MODEL, cwd="/tmp"):
        """
        Get the model version.
        @keyword model_settings: cancer model settings
        @keyword cwd: working directory
        """
        try:
            process = Popen(
                [os.path.join(model['HOME'], model['EXE']), "-v"],
                cwd=cwd,
                stdout=PIPE,
                stderr=PIPE,
                env=settings.FORTRAN_ENV)

            (outs, errs) = process.communicate(timeout=settings.FORTRAN_TIMEOUT)   # timeout in seconds
            exit_code = process.wait()

            if exit_code == 0:
                return outs.decode("utf-8").replace('.exe', '').replace('\n', '')
            else:
                logger.error(outs)
                errs = errs.decode("utf-8").replace('\n', '')
                logger.error(errs)
                raise ModelError(errs)
        except TimeoutExpired as to:
            process.terminate()
            logger.error(model.get('NAME', "")+" PROCESS TIMED OUT.")
            logger.error(to)
            raise TimeOutException()
        except Exception as e:
            logger.error(model.get('NAME', "")+' PROCESS EXCEPTION: '+cwd)
            logger.error(e)
            raise

    @classmethod
    def run(cls, request, bat_file, model_opts, model_params, 
            param_file=None, cwd="/tmp", niceness=0, name="", model=settings.BC_MODEL):
        """
        Run a process.
        @param request: HTTP request
        @param bat_file: batch file path
        @param model_opts: fortran model options
        @param model_params: fortran model parameters
        @param param_file: settings file name
        @keyword cwd: working directory
        @keyword niceness: niceness value
        @keyword name: log name for calculation, e.g. REMAINING LIFETIME
        """
        cancer_rates = model_params.cancer_rates
        cmd = [os.path.join(model['HOME'], model['EXE'])]
        mname = str(model.get('NAME', ""))
        if param_file is not None:
            cmd.extend(["-s", param_file])
        if mname == "BC":
            cmd.extend(["-e", os.path.join(model["HOME"], 'Data', model_params.ethnicity.get_filename())])

        out = model_opts.out
        cmd.extend(model_opts.get_cmd_line_opts())
        cmd.extend([bat_file, model['INCIDENCE'] + cancer_rates + ".nml"])

        start = time.time()
        try:
            try:
                os.remove(os.path.join(cwd, out))  # ensure output file doesn't exist
            except OSError:
                pass

            # logger.debug(' '.join(cmd))
            process = Popen(
                cmd,
                cwd=cwd,
                stdout=PIPE,
                stderr=PIPE,
                env=settings.FORTRAN_ENV,
                preexec_fn=lambda: os.nice(niceness) and
                resource.setrlimit(resource.RLIMIT_STACK, (resource.RLIM_INFINITY, resource.RLIM_INFINITY)))

            (outs, errs) = process.communicate(timeout=settings.FORTRAN_TIMEOUT)   # timeout in seconds
            exit_code = process.wait()

            if exit_code == 0:
                with open(os.path.join(cwd, out), 'r') as result_file:
                    data = result_file.read()
                logger.info(
                    f"{mname} {name} CALCULATION: user={request.user.id}; "
                    f"elapsed time={time.time() - start}")
                return data
            else:
                logger.error(f"EXIT CODE ({out.replace('can_', '')}): {exit_code}")
                logger.error(outs)
                errs = errs.decode("utf-8").replace('\n', '')
                logger.error(errs)
                raise ModelError(errs)
        except TimeoutExpired as to:
            process.terminate()
            logger.error(f"{mname} PROCESS TIMED OUT.")
            logger.error(to)
            raise TimeOutException()
        except Exception as e:
            logger.error(f"{mname} PROCESS EXCEPTION: {cwd}")
            logger.error(e)
            raise
