""" Risk and mutation probability calculations, for details
see https://github.com/CCGE-BOADICEA/boadicea/wiki/Cancer-Risk-Calculations"""
from bws import pedigree
from bws.cancer import Cancer, Cancers, CanRiskGeneticTests, BWSGeneticTests
from bws.pedigree import Male, Female, BwaPedigree, CanRiskPedigree
from bws.exceptions import TimeOutException, ModelError
from collections import OrderedDict
from copy import deepcopy
from django.conf import settings
from django.contrib.auth.models import User, AnonymousUser
from django.http.request import HttpRequest
from rest_framework.request import Request
from subprocess import Popen, PIPE, TimeoutExpired
import logging
import os
import resource
import tempfile
import time
from rest_framework.exceptions import ValidationError


logger = logging.getLogger(__name__)


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

    def _get_prs(self):
        """
        Get the prs.
        @return: prs
        """
        return self.predictions.prs

    def _get_name(self):
        return ""

    def _type(self):
        """ Returns the type of risk as the class name. """
        return self.__class__.__name__

    def get_risk(self):
        """
        Calulate the risk and return the parsed output as a list.
        @return: list of risks for each age
        """
        pedi = self._get_pedi()
        ped_file = pedi.write_pedigree_file(file_type=pedigree.CANCER_RISKS,
                                            risk_factor_code=self._get_risk_factor_code(),
                                            prs=self._get_prs(),
                                            filepath=os.path.join(self.predictions.cwd, self._type()+"_risk.ped"),
                                            model_settings=self.predictions.model_settings)
        bat_file = pedi.write_batch_file(pedigree.CANCER_RISKS, ped_file,
                                         filepath=os.path.join(self.predictions.cwd, self._type()+"_risk.bat"),
                                         mutation_freq=self.predictions.mutation_frequency,
                                         sensitivity=self.predictions.mutation_sensitivity,
                                         calc_ages=self.risk_age)
        risks = Predictions.run(self.predictions.request, pedigree.CANCER_RISKS, bat_file,
                                cancer_rates=self.predictions.cancer_rates, cwd=self.predictions.cwd,
                                niceness=self.predictions.niceness, name=self._get_name(),
                                model=self.predictions.model_settings)
        return self._parse_risks_output(risks)

    def _parse_risks_output(self, risks):
        """
        Parse computed cancer risk results.
        @param risks: cancer risks text from fortran output
        @return: list of containing dictionaries of the risk results for each age
        """
        lines = risks.split(sep="\n")
        risks_arr = []
        version = None
        for line in lines:
            if pedigree.BLANK_LINE.match(line):
                continue
            if line.startswith('#Version:'):
                version = line.replace('#Version:', '').strip()
            elif not line.startswith('#'):
                parts = line.split(sep=",")
                if self.predictions.model_settings['NAME'] == 'BC':
                    risks_arr.append(OrderedDict([
                        ("age", int(parts[0])),
                        ("breast cancer risk", {
                            "decimal": float(parts[1]),
                            "percent": float(parts[2])
                        }),
                        ("ovarian cancer risk", {
                            "decimal": float(parts[3]),
                            "percent": float(parts[4])
                        })
                    ]))
                else:
                    risks_arr.append(OrderedDict([
                        ("age", int(parts[0])),
                        ("ovarian cancer risk", {
                            "decimal": float(parts[1]),
                            "percent": float(parts[2])
                        })
                    ]))
        return risks_arr, version


class RemainingLifetimeRisk(Risk):

    def _get_name(self):
        return "REMAINING LIFETIME"


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

        if t.sex() is "M":
            new_t = Male(t.famid, t.name, t.pid, "", "", target=t.target, dead="0",
                         age=t.age, yob=t.yob, cancers=cancers, gtests=gtests)
        else:
            new_t = Female(t.famid, t.name, t.pid, "", "", target=t.target, dead="0",
                           age=t.age, yob=t.yob, cancers=cancers, gtests=gtests)

        if self.predictions.model_settings['NAME'] == 'BC':
            return BwaPedigree(people=[new_t])
        else:
            return CanRiskPedigree(people=[new_t])

    def _get_risk_factor_code(self):
        return '0'

    def _get_prs(self):
        return None

    def _get_name(self):
        return "REMAINING LIFETIME BASELINE"


class RangeRisk(Risk):
    """
    Calculate risk over a time range, e.g. lifetime or 10 year.
    These should only be calculated for unaffected probands.
    All genetic and epidemiological risk factors should remain unchanged.
    Other family members remain unchanged. (AA email 16/6/2017)
    """
    def __init__(self, predictions, current_age, risk_age, name):
        """
        Run cancer risk and mutation probability prediction calculations.
        @param predictions: L{Predictions} used in prediction calculations
        @param current_age: current age of the proband
        @param risk_age: age to which risk is calculated
        @param name: name of the risk calculation, e.g. LIFETIME
        """
        super().__init__(predictions)
        self.current_age = current_age
        self.risk_age = risk_age
        self.name = name

    def get_risk(self):
        t = self.predictions.pedi.get_target()
        if t.cancers.is_cancer_diagnosed():   # not calculated for affected indivual's
            return (None, None)
        return super().get_risk()

    def _get_pedi(self):
        new_pedi = deepcopy(self.predictions.pedi)
        t = new_pedi.get_target()
        t.age = self.current_age
        return new_pedi

    def _get_name(self):
        return self.name


class RangeRiskBaseline(RangeRisk):
    """
    Calculate baseline risk over a time range.
    """
    def _get_pedi(self):
        t = super()._get_pedi().get_target()
        if t.sex() is "M":
            new_t = Male(t.famid, t.name, t.pid, "", "", target=t.target,
                         dead="0", age=t.age, yob=t.yob, cancers=t.cancers,
                         gtests=t.gtests)
        else:
            new_t = Female(t.famid, t.name, t.pid, "", "", target=t.target,
                           dead="0", age=t.age, yob=t.yob, cancers=t.cancers,
                           gtests=t.gtests)

        if self.predictions.model_settings['NAME'] == 'BC':
            return BwaPedigree(people=[new_t])
        else:
            return CanRiskPedigree(people=[new_t])

    def _get_risk_factor_code(self):
        return '0'

    def _get_prs(self):
        return None


class Predictions(object):

    def __init__(self, pedi,
                 mutation_frequency=settings.BC_MODEL['MUTATION_FREQUENCIES']["UK"],
                 mutation_sensitivity=settings.BC_MODEL['GENETIC_TEST_SENSITIVITY'],
                 cancer_rates=settings.BC_MODEL['CANCER_RATES'].get("UK"),
                 risk_factor_code=0, prs=None, cwd=None, request=Request(HttpRequest()),
                 run_risks=True, model_settings=settings.BC_MODEL, calcs=None):
        """
        Run cancer risk and mutation probability prediction calculations.
        @param pedi: L{Pedigree} used in prediction calculations
        @keyword mutation_frequrency: mutation frequencies used in model
        @keyword mutation_sensitivity: mutation search sensitivities
        @keyword cancer_rates: cancer incidence rates used in risk calculation
        @keyword risk_factor_code: risk factor code
        @keyword prs: polygenic risk alpha & beta values calculated from VCF file
        @keyword cwd: working directory
        @keyword request: HTTP request
        @keyword run_risks: run risk calculations, default True
        @keyword model_settings: cancer model settings
        @keyword calcs: list of calculations to run, e.g. ['carrier_probs', 'remaining_lifetime']
        """
        self.pedi = pedi
        self.mutation_frequency = mutation_frequency
        self.mutation_sensitivity = mutation_sensitivity
        self.cancer_rates = cancer_rates
        self.request = request
        self.cwd = cwd
        self.risk_factor_code = risk_factor_code
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
        self.niceness = Predictions._get_niceness(self.pedi)
        start = time.time()

        # mutation probability calculation
        if self.pedi.is_carrier_probs_viable() and self.is_calculate('carrier_probs'):
            ped_file = self.pedi.write_pedigree_file(file_type=pedigree.MUTATION_PROBS,
                                                     risk_factor_code=self.risk_factor_code,
                                                     prs=self.prs,
                                                     filepath=os.path.join(self.cwd, "test_prob.ped"),
                                                     model_settings=self.model_settings)
            bat_file = self.pedi.write_batch_file(pedigree.MUTATION_PROBS, ped_file,
                                                  filepath=os.path.join(self.cwd, "test_prob.bat"),
                                                  mutation_freq=self.mutation_frequency,
                                                  sensitivity=self.mutation_sensitivity)
            probs = self.run(self.request, pedigree.MUTATION_PROBS, bat_file, cancer_rates=self.cancer_rates,
                             cwd=self.cwd, niceness=self.niceness, model=self.model_settings)
            self.mutation_probabilties, self.version = self._parse_probs_output(probs, self.model_settings)

        # cancer risk calculation
        if self.pedi.is_risks_calc_viable():
            # remaining lifetime risk
            if self.is_calculate("remaining_lifetime"):
                self.cancer_risks, self.version = RemainingLifetimeRisk(self).get_risk()
                self.baseline_cancer_risks, _v = RemainingLifetimeBaselineRisk(self).get_risk()

            # lifetime risk
            if self.is_calculate("lifetime"):
                self.lifetime_cancer_risk, _v = RangeRisk(self, 20, 80, "LIFETIME").get_risk()
                if self.lifetime_cancer_risk is not None:
                    self.baseline_lifetime_cancer_risk, _v = RangeRiskBaseline(self, 20, 80,
                                                                               "LIFETIME BASELINE").get_risk()

            # ten year risk
            if self.is_calculate("ten_year"):
                self.ten_yr_cancer_risk, _v = RangeRisk(self, 40, 50, "10 YR RANGE").get_risk()
                if self.ten_yr_cancer_risk is not None:
                    self.baseline_ten_yr_cancer_risk, _v = RangeRiskBaseline(self, 40, 50,
                                                                             "10YR RANGE BASELINE").get_risk()

        if not isinstance(self.request.user, AnonymousUser):
            u = User.objects.get(username=self.request.user)
            job_title = u.userdetails.job_title
            country = u.userdetails.country.name
        else:
            job_title = 'AnonymousUser'
            country = 'unknown'
        logger.info(self.model_settings.get('NAME', "") + " CALCULATIONS: " +
                    "job=" + job_title +
                    "; country=" + country +
                    "; elapsed time=" + str(time.time() - start) +
                    "; pedigree size=" + str(len(self.pedi.people)) +
                    "; version=" + (getattr(self, "version", "N/A")))

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
        if len(siblings) > 0:
            return 19

        pedigree_size = len(pedi.people)
        niceness = int(pedigree_size/factor)
        if niceness > 19:
            niceness = 19
        return niceness

    @classmethod
    def run(cls, request, process_type, bat_file, cancer_rates="UK", cwd="/tmp", niceness=0, name="",
            model=settings.BC_MODEL):
        """
        Run a process.
        @param request: HTTP request
        @param process_type: either pedigree.MUTATION_PROBS or pedigree.CANCER_RISKS.
        @param bat_file: batch file path
        @keyword cancer_rates: cancer incidence rates used in risk calculation
        @keyword cwd: working directory
        @keyword niceness: niceness value
        @keyword name: log name for calculation, e.g. REMAINING LIFETIME
        """
        if process_type == pedigree.MUTATION_PROBS:
            prog = os.path.join(model['HOME'], model['PROBS_EXE'])
            out = "can_probs"
        else:
            prog = os.path.join(model['HOME'], model['RISKS_EXE'])
            out = "can_risks"

        start = time.time()
        try:
            try:
                os.remove(os.path.join(cwd, out+".out"))  # ensure output file doesn't exist
            except OSError:
                pass

            process = Popen(
                [prog,
                 bat_file,  # "Sample_Pedigrees/risks_single_person.bat",
                 os.path.join(model['HOME'], "Data/incidence_rates_" + cancer_rates + ".nml"),
                 out+".out"],
                cwd=cwd,
                stdout=PIPE,
                stderr=PIPE,
                env=settings.FORTRAN_ENV,
                preexec_fn=lambda: os.nice(niceness) and
                resource.setrlimit(resource.RLIMIT_STACK, (resource.RLIM_INFINITY, resource.RLIM_INFINITY)))

            (outs, errs) = process.communicate(timeout=settings.FORTRAN_TIMEOUT)   # timeout in seconds
            exit_code = process.wait()

            if exit_code == 0:
                with open(os.path.join(cwd, out+".out"), 'r') as result_file:
                    data = result_file.read()
                logger.info(model.get('NAME', "") + " " +
                            ("MUTATION PROBABILITY" if process_type == pedigree.MUTATION_PROBS else "RISK ") +
                            name + " CALCULATION: " + "elapsed time=" + str(time.time() - start))
                return data
            else:
                logger.error("EXIT CODE ("+out.replace('can_', '')+"): "+str(exit_code))
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

    def _parse_probs_output(self, probs, model_settings):
        """
        Parse computed mutation carrier probability results.
        @param probs: mutation probaility text from fortran output
        @param model_settings: cancer model settings
        @return: list of containing dictionaries of the mutaion probability results
        """
        probs_arr = []
        version = None
        for line in probs.splitlines():
            if line.startswith('#Version:'):
                version = line.replace('#Version:', '').strip()
            elif not line.startswith('#'):
                parts = line.strip().split(sep=",")
                probs_arr.append({"no mutation": {"decimal": float(parts[0]), "percent": float(parts[1])}})
                for i, gene in enumerate(model_settings['GENES']):
                    probs_arr.append({gene:
                                      {"decimal": float(parts[((i*2)+2)]),
                                       "percent": float(parts[(i*2)+3])}})
        return probs_arr, version
