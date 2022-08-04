""" Risk and mutation probability calculations, for details
see https://github.com/CCGE-BOADICEA/boadicea/wiki/Cancer-Risk-Calculations"""
from collections import OrderedDict
import logging
import os
import resource
from subprocess import Popen, PIPE, TimeoutExpired
import tempfile
import time

from django.conf import settings
from django.http.request import HttpRequest
from rest_framework.exceptions import ValidationError, NotAcceptable
from rest_framework.request import Request

from bws import pedigree
from bws.cancer import Cancer, Cancers, CanRiskGeneticTests, BWSGeneticTests
from bws.exceptions import TimeOutException, ModelError
from bws.pedigree import Male, Female, BwaPedigree, CanRiskPedigree
import re


logger = logging.getLogger(__name__)

REGEX_ALPHANUM_COMMAS = re.compile("^([\\w,]+)$")


class ModelOpts():

    """
    Fortran cancer model options
    -o , --output= write results to a file [defaults to the stdout].
    -p, --probs calculates pathogenic carrier probabilities.
    -rj, --riskJ 10-year cancer risk calculation (NHS protocol for young women at high risk)
    -rl, --riskL lifetime cancer risk calculation (proband's age set at 20y, censor age set at 80y)
    -rr, --riskR cancer risk calculation is performed (proband's age and censor age set in input files)
    -ry, --riskY 10-year cancer risk calculation (proband's age set at 40y, censor age set at 50y)
    -s , --settings= use values of the model parameters in the specified file [default model settings if absent].
    -v, --version print the version.
    """
    def __init__(self, out="predictions.txt", probs=True, rj=True, rl=True, rr=True, ry=True):
        self.out = out
        self.probs = probs
        self.rj = rj
        self.rl = rl
        self.rr = rr
        self.ry = ry

    def get_cmd_line_opts(self):
        cmd = ["-o", self.out]
        if self.probs:
            cmd.extend(["-p"])
        if self.rj:
            cmd.extend(["-rj"])
        if self.rl:
            cmd.extend(["-rl"])
        if self.rr:
            cmd.extend(["-rr"])
        if self.ry:
            cmd.extend(["-ry"])
        return cmd

    @classmethod
    def factory(cls, calc):
        """
        Generate ModelOpts from web-service validated data
        @keyword calc: Predictions
        @return: ModelOpts
        """
        t = calc.pedi.get_target()
        is_alive = (t.dead != "1")
        is_cancer_diagnosed = t.cancers.is_cancer_diagnosed()
        is_risks_calc_viable = calc.pedi.is_risks_calc_viable() and is_alive and t.sex() == "F"
        is_carr_probs_viable = calc.pedi.is_carrier_probs_viable()
        is_nhs_tenyr_reqd = (int(t.age) < 50 and calc.model_settings['NAME'] == 'BC')

        mname = str(calc.model_settings.get('NAME', ""))
        return ModelOpts(out=mname+"_predictions.txt",
                         probs=(is_carr_probs_viable and calc.is_calculate('carrier_probs')),
                         rj=is_risks_calc_viable and is_nhs_tenyr_reqd,
                         rl=(is_risks_calc_viable and calc.is_calculate("lifetime") and not is_cancer_diagnosed),
                         rr=is_risks_calc_viable,
                         ry=(is_risks_calc_viable and calc.is_calculate("ten_year") and not is_cancer_diagnosed))


class ModelParams():

    def __init__(self, population="UK", mutation_frequency=settings.BC_MODEL['MUTATION_FREQUENCIES']["UK"],
                 mutation_sensitivity=settings.BC_MODEL['GENETIC_TEST_SENSITIVITY'],
                 cancer_rates=settings.BC_MODEL['CANCER_RATES'].get("UK")):
        """
        Cancer risk model parameters and population.
        @keyword population: population setting
        @keyword mutation_frequrency: mutation frequencies used in model
        @keyword mutation_sensitivity: mutation search sensitivities
        @keyword cancer_rates: cancer incidence rates used in risk calculation
        """
        self.population = population        # e.g. Ashkenazi, UK, Iceland
        self.cancer_rates = cancer_rates
        self.mutation_frequency = mutation_frequency
        self.mutation_sensitivity = mutation_sensitivity
        self.isashk = settings.REGEX_ASHKN.match(population)

    @classmethod
    def factory(cls, data, model_settings):
        """
        Generate ModelParams from web-service validated data
        @keyword data: validated request data
        @keyword model_settings: model settings
        @return: ModelParams
        """
        population = data.get('mut_freq', 'UK')
        crates = model_settings['CANCER_RATES'].get(data.get('cancer_rates'))

        if population != 'Custom':
            mut_freq = model_settings['MUTATION_FREQUENCIES'][population]
        else:
            mut_freq = {}
            for gene in model_settings['GENES']:
                try:
                    mut_freq[gene] = float(data.get(gene.lower() + '_mut_frequency'))
                except TypeError:
                    raise NotAcceptable("Invalid mutation frequency for " + gene + ".")

        gts = model_settings['GENETIC_TEST_SENSITIVITY']
        mut_sens = {
            k: float(data.get(k.lower() + "_mut_sensitivity", gts[k]))
            for k in gts.keys()
        }
        return ModelParams(population, cancer_rates=crates, mutation_frequency=mut_freq, mutation_sensitivity=mut_sens)


class Risk(object):

    def __init__(self, predictions):
        """
        Run cancer risk and mutation probability prediction calculations.
        @param predictions: L{Predictions} used in prediction calculations
        """
        assert isinstance(predictions, Predictions), "%r is not a Predictions" % predictions
        self.predictions = predictions
        self.risk_age = None    # if none uses the default ages to calculate risk at

    def _get_pedi(self):
        """
        Get the pedigree.
        @return: L{Pedigree}
        """
        return self.predictions.pedi

    def _get_risk_factor_code(self):
        """
        Get the risk factor code.
        @return: risk factor code
        """
        return self.predictions.risk_factor_code

    def _get_hgt(self):
        """
        Get the height.
        @return: height
        """
        return self.predictions.hgt

    def _get_mutation_frequency(self):
        """
        Get the mutation frequencies.
        @return: mutation frequencies
        """
        return self.predictions.model_params.mutation_frequency

    def _get_prs(self):
        """
        Get the prs.
        @return: prs
        """
        return self.predictions.prs

    def _get_name(self):
        return "RISK AND MUTATION CARRIER PREDICTIONS"

    def _type(self):
        """ Returns the type of risk as the class name. """
        return self.__class__.__name__

    def get_risk(self, model_opts):
        """
        Calculate the risk and return the parsed output as a list.
        @return: list of risks for each age
        """
        pedi = self._get_pedi()
        pred = self.predictions
        ped_file = pedi.write_pedigree_file(file_type=pedigree.CANCER_RISKS,
                                            risk_factor_code=self._get_risk_factor_code(),
                                            hgt=self._get_hgt(),
                                            prs=self._get_prs(),
                                            filepath=os.path.join(pred.cwd, self._type()+"_risk.ped"),
                                            model_settings=pred.model_settings)
        bat_file = pedi.write_batch_file(pedigree.CANCER_RISKS, ped_file,
                                         filepath=os.path.join(pred.cwd, self._type()+"_risk.bat"),
                                         model_settings=pred.model_settings,
                                         calc_ages=self.risk_age)
        params = pedi.write_param_file(filepath=os.path.join(pred.cwd, self._type()+"_risk.params"),
                                       model_settings=pred.model_settings,
                                       mutation_freq=self._get_mutation_frequency(),
                                       isashk=pred.model_params.isashk,
                                       sensitivity=pred.model_params.mutation_sensitivity)
        risks = Predictions.run(self.predictions.request, pedigree.CANCER_RISKS, bat_file,
                                model_opts=model_opts,
                                params=params,
                                cancer_rates=pred.model_params.cancer_rates, cwd=pred.cwd,
                                niceness=pred.niceness, name=self._get_name(),
                                model=pred.model_settings)
        return self._parse_risks_output(risks, model_opts)

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
        model_settings = self.predictions.model_settings
        ctype = "breast" if model_settings['NAME'] == 'BC' else "ovarian"

        for _idx, line in enumerate(lines):
            if line.startswith('##'):
                rr, rl, ry, rj, mp = False, False, False, False, False
            elif 'Age' in line or pedigree.BLANK_LINE.match(line):
                continue

            if rr:
                parts = line.split(sep=",")
                rr_arr.append(OrderedDict([
                    ("age", int(parts[0])),
                    (ctype+" cancer risk", {
                        "decimal": float(parts[1]),
                        "percent": round(float(parts[1])*100, 1)
                    })
                ]))
            elif rl:
                parts = line.split(sep=",")
                rl_arr.append(OrderedDict([
                    ("age", int(parts[1])),
                    (ctype+" cancer risk", {
                        "decimal": float(parts[2]),
                        "percent": round(float(parts[2])*100, 1)
                    })
                ]))
            elif ry and ry_arr is not None:
                parts = line.split(sep=",")
                ry_arr.append(OrderedDict([
                    ("age", int(parts[1])),
                    (ctype+" cancer risk", {
                        "decimal": float(parts[2]),
                        "percent": round(float(parts[2])*100, 1)
                    })
                ]))
            elif rj:
                parts = line.split(sep=",")
                rj_arr.append(OrderedDict([
                    ("age", int(parts[1])),
                    (ctype+" cancer risk", {
                        "decimal": float(parts[2]),
                        "percent": round(float(parts[2])*100, 1)
                    })
                ]))
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

        mp_arr = self._parse_probs_output(mp_lines, model_settings) if model_opts.probs else None
        return rl_arr, rr_arr, ry_arr, rj_arr, mp_arr

    def _parse_probs_output(self, probs, model_settings):
        """
        Parse computed mutation carrier probability results.
        @param probs: mutation probaility text from fortran output
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


class RemainingLifetimeBaselineRisk(Risk):
    """
    Get the baseline risks: the purpose of the baseline risk is to show the risk to
    an equivalent random woman from the population without any information on risk of
    genetic factors (i.e. based on population incidences only). The only aspects
    observed are her age (she needs to alive..!) year of birth and age at cancer diagnosis
    if affected. (ACA: 16/2/2017 email)
    """
    def _get_pedi(self):
        t = self.predictions.pedi.get_target()
        if t.cancers.is_cancer_diagnosed():
            cancers = Cancers(bc1=Cancer(t.cancers.diagnoses.bc1.age), bc2=Cancer(), oc=Cancer(),
                              prc=Cancer(), pac=Cancer())
        else:
            cancers = Cancers()

        if self.predictions.model_settings['NAME'] == 'BC':
            gtests = BWSGeneticTests.default_factory()
        else:
            gtests = CanRiskGeneticTests.default_factory()

        if t.sex() == "M":
            new_t = Male(t.famid, t.name, t.pid, "", "", target=t.target, dead=t.dead,
                         age=t.age, yob=t.yob, cancers=cancers, gtests=gtests)
        else:
            new_t = Female(t.famid, t.name, t.pid, "", "", target=t.target, dead=t.dead,
                           age=t.age, yob=t.yob, cancers=cancers, gtests=gtests)

        if self.predictions.model_settings['NAME'] == 'BC':
            return BwaPedigree(people=[new_t])
        else:
            return CanRiskPedigree(people=[new_t])

    def _get_risk_factor_code(self):
        return '0'

    def _get_hgt(self):
        return -1

    def _get_prs(self):
        return None

    def _get_name(self):
        return "REMAINING LIFETIME BASELINE"


class RiskBaseline(Risk):
    """
    Calculate baseline risk over a time range.
    """
    def _get_pedi(self):
        t = super()._get_pedi().get_target()
        cancers = Cancers()
        if self.predictions.model_settings['NAME'] == 'BC':
            gtests = BWSGeneticTests.default_factory()
        else:
            gtests = CanRiskGeneticTests.default_factory()

        if t.sex() == "M":
            new_t = Male(t.famid, t.name, t.pid, "", "", target=t.target,
                         dead=t.dead, age=t.age, yob=t.yob, cancers=cancers,
                         gtests=gtests)
        else:
            new_t = Female(t.famid, t.name, t.pid, "", "", target=t.target,
                           dead=t.dead, age=t.age, yob=t.yob, cancers=cancers,
                           gtests=gtests)

        if self.predictions.model_settings['NAME'] == 'BC':
            return BwaPedigree(people=[new_t])
        else:
            return CanRiskPedigree(people=[new_t])

    def _get_risk_factor_code(self):
        return '0'

    def _get_hgt(self):
        return -1

    def _get_prs(self):
        return None

    def _get_name(self):
        return "BASELINE RISK PREDICTIONS"


class Predictions(object):

    def __init__(self, pedi, model_params=ModelParams(),
                 risk_factor_code=0, hgt=-1, prs=None, cwd=None, request=Request(HttpRequest()),
                 run_risks=True, model_settings=settings.BC_MODEL, calcs=None):
        """
        Run cancer risk and mutation probability prediction calculations.
        @param pedi: L{Pedigree} used in prediction calculations
        @keyword model_params: model parameters
        @keyword risk_factor_code: risk factor code
        @keyword prs: polygenic risk alpha & beta values calculated from VCF file
        @keyword cwd: working directory
        @keyword request: HTTP request
        @keyword run_risks: run risk calculations, default True
        @keyword model_settings: cancer model settings
        @keyword calcs: list of calculations to run, e.g. ['carrier_probs', 'remaining_lifetime']
        """
        self.pedi = pedi
        self.model_params = model_params
        self.request = request
        self.cwd = cwd
        self.risk_factor_code = risk_factor_code
        self.hgt = hgt
        self.prs = prs
        self.model_settings = model_settings
        self.calcs = self.model_settings['CALCS'] if calcs is None else calcs

        # check calculations are in the allowed list of calculations
        for c in self.calcs:
            if c not in settings.ALLOWED_CALCS:
                raise ValidationError("Unknown calculation requested: "+c)

        if cwd is None:
            self.cwd = tempfile.mkdtemp(prefix=str(request.user)+"_", dir="/tmp")
        if isinstance(risk_factor_code, int):
            self.risk_factor_code = str(risk_factor_code)
        if run_risks:
            self.run_risks()

    def is_calculate(self, calc):
        '''
        Determine if a calculation is to be run.
        @param calc: calculation name, e.g. 'carrier_probs', 'remaining_lifetime'
        '''
        return True if len(self.calcs) == 0 else (calc in self.calcs)

    def run_risks(self):
        ''' Run risk and mutation probability calculations '''
        self.version = Predictions.get_version(model=self.model_settings, cwd=self.cwd)
        self.niceness = Predictions._get_niceness(self.pedi)
        start = time.time()
        model_opts = ModelOpts.factory(self)

        rl, rr, ry, rj, mp = Risk(self).get_risk(model_opts)
        if rl is not None:
            self.lifetime_cancer_risk = rl
        if rr is not None:
            self.cancer_risks = rr
        if ry is not None:
            self.ten_yr_cancer_risk = ry
        if rj is not None:
            self.ten_yr_nhs_protocol = rj
        if mp is not None:
            self.mutation_probabilties = mp

        # remaining lifetime baseline
        if rr is not None:
            _rl, rr, _ry, _rj, _mp = RemainingLifetimeBaselineRisk(self).get_risk(ModelOpts(probs=False,
                                                                                            rj=False, rl=False,
                                                                                            rr=True, ry=False))
            if rr is not None:
                self.baseline_cancer_risks = rr

            # baseline lifetime cancer risk and baseline 10-year cancer risk
            rl, _rr, ry, _rj, _mp = RiskBaseline(self).get_risk(ModelOpts(probs=False, rj=False, rl=True,
                                                                          rr=False, ry=True))
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
    def get_version(cls, model=settings.BC_MODEL, cwd="/tmp"):
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
    def run(cls, request, process_type, bat_file, model_opts, params=None, cancer_rates="UK", cwd="/tmp",
            niceness=0, name="", model=settings.BC_MODEL):
        """
        Run a process.
        @param request: HTTP request
        @param process_type: either pedigree.MUTATION_PROBS or pedigree.CANCER_RISKS.
        @param bat_file: batch file path
        @param model_opts: fortran model options
        @param params: model parameters
        @keyword cancer_rates: cancer incidence rates used in risk calculation
        @keyword cwd: working directory
        @keyword niceness: niceness value
        @keyword name: log name for calculation, e.g. REMAINING LIFETIME
        """
        cmd = [os.path.join(model['HOME'], model['EXE'])]

        if params is not None:
            cmd.extend(["-s", params])

        out = model_opts.out
        cmd.extend(model_opts.get_cmd_line_opts())
        cmd.extend([bat_file, model['INCIDENCE'] + cancer_rates + ".nml"])
        mname = str(model.get('NAME', ""))

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
