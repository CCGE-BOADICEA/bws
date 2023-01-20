'''
Model options and parameters

© 2022 Cambridge University
SPDX-FileCopyrightText: 2022 Cambridge University
SPDX-License-Identifier: GPL-3.0-or-later
'''
from django.conf import settings
from bws.risk_factors.ethnicity import UKBioBankEthnicty


class ModelOpts():

    """
    Defines the Fortran cancer risk model options
    -o , --output= write results to a file [defaults to the stdout].
    -p, --probs calculates pathogenic carrier probabilities.
    -rj, --riskJ 10-year cancer risk calculation (NHS protocol for young women at high risk)
    -rl, --riskL lifetime cancer risk calculation (proband's age set at 20y, censor age set at 80y)
    -rr, --riskR cancer risk calculation is performed (proband's age and censor age set in input files)
    -ry, --riskY 10-year cancer risk calculation (proband's age set at 40y, censor age set at 50y)
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
        mname = calc.model_settings.get('NAME', "")
        t = calc.pedi.get_target()
        is_alive = (t.dead != "1")
        is_cancer_diagnosed = t.cancers.is_cancer_diagnosed()
        is_risks_calc_viable = calc.pedi.is_risks_calc_viable() and is_alive and t.sex() == "F"
        is_carr_probs_viable = calc.pedi.is_carrier_probs_viable()
        return ModelOpts(out=mname+"_predictions.txt",
                         probs=(is_carr_probs_viable and calc.is_calculate('carrier_probs')),
                         rj=(is_risks_calc_viable and (mname == 'BC' and int(t.age) < 50) and not is_cancer_diagnosed),
                         rl=(is_risks_calc_viable and calc.is_calculate("lifetime") and not is_cancer_diagnosed),
                         rr=is_risks_calc_viable,
                         ry=(is_risks_calc_viable and calc.is_calculate("ten_year") and not is_cancer_diagnosed))


class ModelParams():

    def __init__(self, population="UK", mutation_frequency=settings.BC_MODEL['MUTATION_FREQUENCIES']["UK"],
                 mutation_sensitivity=settings.BC_MODEL['GENETIC_TEST_SENSITIVITY'],
                 cancer_rates=settings.BC_MODEL['CANCER_RATES'].get("UK"),
                 ethnicity=UKBioBankEthnicty()):
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
        self.ethnicity = ethnicity

    @classmethod
    def factory(cls, data, model_settings):
        """
        Generate ModelParams from web-service validated data
        @keyword data: validated request data
        @keyword model_settings: model settings
        @return: ModelParams
        """
        population = data.get('mut_freq', 'UK')
        gts = model_settings['GENETIC_TEST_SENSITIVITY']
        mut_sens = {k: float(data.get(k.lower() + "_mut_sensitivity", gts[k])) for k in gts.keys()}
        return ModelParams(population,
                           cancer_rates=model_settings['CANCER_RATES'].get(data.get('cancer_rates')),
                           mutation_frequency=model_settings['MUTATION_FREQUENCIES'][population],
                           mutation_sensitivity=mut_sens)
