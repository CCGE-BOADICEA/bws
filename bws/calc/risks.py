"""
Risk calculations

© 2023 University of Cambridge
SPDX-FileCopyrightText: 2023 University of Cambridge
SPDX-License-Identifier: GPL-3.0-or-later
"""
from bws.cancer import Cancer, Cancers, CanRiskGeneticTests, BWSGeneticTests
from bws.pedigree import Male, Female, BwaPedigree, CanRiskPedigree


class Risk(object):

    def __init__(self, predictions):
        """
        Run cancer risk and mutation probability prediction calculations.
        @param predictions: L{Predictions} used in prediction calculations
        """
        self.predictions = predictions
        self.risk_age = None    # if none uses the default ages to calculate risk at

    def get_pedigree(self):
        """
        Get the pedigree.
        @return: L{Pedigree}
        """
        return self.predictions.pedi

    def get_risk_factor_code(self):
        """
        Get the risk factor code.
        @return: risk factor code
        """
        return self.predictions.risk_factor_code

    def get_hgt(self):
        """
        Get the height.
        @return: height
        """
        return self.predictions.hgt

    def get_md(self):
        """
        Get the mammographic density 
        @return: mammographic density
        """
        return self.predictions.mdensity

    def get_mutation_frequency(self):
        """
        Get the mutation frequencies.
        @return: mutation frequencies
        """
        return self.predictions.model_params.mutation_frequency

    def get_prs(self):
        """
        Get the prs.
        @return: prs
        """
        return self.predictions.prs

    def get_name(self):
        return "RISK AND MUTATION CARRIER PREDICTIONS"

    def type(self):
        """ Returns the type of risk as the class name. """
        return self.__class__.__name__


class RemainingLifetimeBaselineRisk(Risk):
    """
    Get the baseline risks: the purpose of the baseline risk is to show the risk to
    an equivalent random woman from the population without any information on risk of
    genetic factors (i.e. based on population incidences only). The only aspects
    observed are her age (she needs to alive..!) year of birth and age at cancer diagnosis
    if affected. (ACA: 16/2/2017 email)
    """
    def get_pedigree(self):
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

    def get_risk_factor_code(self):
        return '0'

    def get_hgt(self):
        return -1

    def get_md(self):
        """
        Get the mammographic density.
        @return: mammographic density
        """
        return None

    def get_prs(self):
        return None

    def get_name(self):
        return "REMAINING LIFETIME BASELINE"


class RiskBaseline(Risk):
    """
    Calculate baseline risk over a time range.
    """
    def get_pedigree(self):
        t = super().get_pedigree().get_target()
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

    def get_risk_factor_code(self):
        return '0'

    def get_hgt(self):
        return -1

    def get_md(self):
        """
        Get the mammographic density.
        @return: mammographic density
        """
        return None

    def get_prs(self):
        return None

    def get_name(self):
        return "BASELINE RISK PREDICTIONS"
