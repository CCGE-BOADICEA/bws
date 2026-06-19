"""
Cancer risk and pathogenic variant probability calculation testing.

© 2023 University of Cambridge
SPDX-FileCopyrightText: 2023 University of Cambridge
SPDX-License-Identifier: GPL-3.0-or-later
"""
import pytest
from copy import deepcopy
import copy
from datetime import date
import filecmp
import os
import shutil
import tempfile

from django.conf import settings
from django.test import TestCase

from bws.calc.calcs import Predictions
from bws.calc.model import ModelParams
from bws.cancer import Cancer, Cancers, CanRiskGeneticTests
from bws.pedigree import Female, BwaPedigree, CanRiskPedigree
from bws.pedigree_file import PedigreeFile


class RiskTests(TestCase):
    """ Calculation testing. """

    def setUp(self):
        ''' Build pedigree data. '''
        self.year = date.today().year

        target = Female("FAM1", "F0", "001", "002", "003", target="1", age="20",
                        yob=str(self.year-20), cancers=Cancers())
        self.pedigree = BwaPedigree(people=[target])
        # parents
        (_father, _mother) = self.pedigree.add_parents(target)

        # canrisk pedigree
        target = Female("FAM1", "F0", "001", "002", "003", target="1", age="20",
                        yob=str(self.year-20), cancers=Cancers(),
                        gtests=CanRiskGeneticTests.default_factory())
        self.canrisk_pedigree = CanRiskPedigree(people=[target])
        (_father, _mother) = self.canrisk_pedigree.add_parents(target, gtests=CanRiskGeneticTests.default_factory())
        _mother.yob = str(self.year-55)
        _mother.age = "55"
        _mother.cancers = Cancers(bc1=Cancer("54"), bc2=Cancer("55"), oc=Cancer(), prc=Cancer(), pac=Cancer())
        self.cwd = tempfile.mkdtemp(prefix="TEST_", dir="/tmp")

    def tearDown(self):
        TestCase.tearDown(self)
        shutil.rmtree(self.cwd)

    @pytest.mark.req_WS_VALIDATION_240
    @pytest.mark.req_WS_VALIDATION_241
    @pytest.mark.req_WS_VALIDATION_242
    def test_risk_calc_viable(self):
        """ Test if the risk calculations are viable. """
        pedigree = deepcopy(self.pedigree)
        target = pedigree.get_target()
        self.assertTrue(pedigree.is_risks_calc_viable())
        self.assertFalse(pedigree.is_risks_calc_viable(allowMale=True))

        # BC1 diagnosis
        target.cancers = Cancers(bc1=Cancer("19"), bc2=Cancer(), oc=Cancer(), prc=Cancer(), pac=Cancer())
        self.assertTrue(pedigree.is_risks_calc_viable())

        # OC diagnosis
        target.cancers = Cancers(bc1=Cancer(), bc2=Cancer(), oc=Cancer("19"), prc=Cancer(), pac=Cancer())
        self.assertFalse(pedigree.is_risks_calc_viable())

    @pytest.mark.req_WS_VALIDATION_243
    @pytest.mark.req_WS_VALIDATION_244
    def test_calculations(self):
        """ Test BC prediction of cancer risk and mutation probability. """
        pedigree = deepcopy(self.pedigree)
        pedigree.validateAll()
        calcs = Predictions(pedigree, cwd=self.cwd)

        # each gene should have a mutation probability plus a result for no mutations
        for mp in calcs.mutation_probabilties:
            key = list(mp.keys())[0]
            self.assertTrue(key in settings.BC_MODEL['GENES'] or key == "no mutation")
        self.assertEqual(len(calcs.mutation_probabilties), len(settings.BC_MODEL['GENES']) + 1)

        # risks calculated at 16 different ages:
        self.assertEqual(len(calcs.cancer_risks), 16)
        self.assertTrue([c.get('age') for c in calcs.cancer_risks] ==
                        [21, 22, 23, 24, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80])

    @pytest.mark.req_WS_VALIDATION_243
    @pytest.mark.req_WS_VALIDATION_244
    def test_calculations2(self):
        """ Test BC prediction of cancer risk and mutation probability with a secondary degree family history. """
        pedigree = deepcopy(self.pedigree)
        target = pedigree.get_target()
        target.age = str(int(target.age) + 43)
        # sister
#         diagnoses = CancerDiagnoses(bc1=Cancer("20"), bc2=Cancer(), oc=Cancer(),
#                                     prc=Cancer(), pac=Cancer())
#         sister = Female("FAM1", "F01", "0011", target.fathid, target.mothid, age="22", yob=str(self.year-23),
#                         cancers=Cancers(diagnoses=diagnoses))
#         pedigree.people.append(sister)

        # parents
        mother = pedigree.get_person(target.mothid)
        mother.yob = str(self.year-84)
        mother.age = "85"
        mother.cancers = Cancers(bc1=Cancer("52"), bc2=Cancer(), oc=Cancer(), prc=Cancer(), pac=Cancer())

        # maternal grandparents
        (_maternal_grandfather, maternal_grandmother) = pedigree.add_parents(mother)

        maternal_grandmother.age = "81"
        maternal_grandmother.yob = "1912"
        maternal_grandmother.dead = "1"
        maternal_grandmother.cancers = Cancers(bc1=Cancer("42"), bc2=Cancer(), oc=Cancer(),
                                               prc=Cancer(), pac=Cancer())

        pedigree.validateAll()
        calcs = Predictions(pedigree, cwd=self.cwd)

        # each gene should have a mutation probability plus a result for no mutations
        for mp in calcs.mutation_probabilties:
            key = list(mp.keys())[0]
            self.assertTrue(key in settings.BC_MODEL['GENES'] or key == "no mutation")
        self.assertEqual(len(calcs.mutation_probabilties), len(settings.BC_MODEL['GENES']) + 1)

        # risks calculated at different ages:
        self.assertEqual(len(calcs.cancer_risks), 9)
        self.assertTrue([c.get('age') for c in calcs.cancer_risks] ==
                        [64, 65, 66, 67, 68, 70, 73, 75, 80])

    @pytest.mark.req_WS_VALIDATION_243
    @pytest.mark.req_WS_VALIDATION_244
    def test_ovarian_calculations(self):
        """ Test OC prediction of cancer risk and mutation probability. """
        target = Female("FAM1", "F0", "001", "002", "003", target="1", age="20",
                        yob=str(self.year-20), cancers=Cancers(),
                        gtests=CanRiskGeneticTests.default_factory())
        pedigree = CanRiskPedigree(people=[target])
        # parents
        (_father, _mother) = pedigree.add_parents(target, gtests=CanRiskGeneticTests.default_factory())
        pedigree.validateAll()
        params = ModelParams(mutation_frequency=settings.OC_MODEL['MUTATION_FREQUENCIES']["UK"],
                             mutation_sensitivity=settings.OC_MODEL['GENETIC_TEST_SENSITIVITY'])
        calcs = Predictions(pedigree, cwd=self.cwd, model_params=params,
                            model_settings=settings.OC_MODEL, calcs=[])

        # each gene should have a mutation probability plus a result for no mutations
        for mp in calcs.mutation_probabilties:
            key = list(mp.keys())[0]
            self.assertTrue(key in settings.OC_MODEL['GENES'] or key == "no mutation")
        self.assertEqual(len(calcs.mutation_probabilties), len(settings.OC_MODEL['GENES']) + 1)

        # risks calculated at 16 different ages:
        self.assertEqual(len(calcs.cancer_risks), 16)
        self.assertTrue([c.get('age') for c in calcs.cancer_risks] ==
                        [21, 22, 23, 24, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80])

    @classmethod
    def run_calc(cls, can, pedigree, calcs=None, cwd=None):
        """ Run rrisk calculations """
        msettings = settings.BC_MODEL if can == 'breast' else settings.OC_MODEL
        mparams = ModelParams(mutation_frequency=msettings['MUTATION_FREQUENCIES']["UK"],
                              mutation_sensitivity=msettings['GENETIC_TEST_SENSITIVITY'])
        calcs = Predictions(pedigree, cwd=cwd, model_settings=msettings, model_params=mparams,
                            calcs=calcs)
        c80 = None
        for c in calcs.cancer_risks:
            if c.get('age') == 80:
                c80 = c
        return c80[can+' cancer risk']['decimal'], calcs.lifetime_cancer_risk[0][can + ' cancer risk']['decimal']

    @pytest.mark.req_WS_VALIDATION_245
    @pytest.mark.req_WS_VALIDATION_246
    def test_lifetime_risk(self):
        """ Test lifetime risk calculation matches the risk to 80 for a 20 year old. """
        # 1. breast cancer risk
        pedigree = deepcopy(self.pedigree)
        t = pedigree.get_target()
        t.age = '20'
        pedigree.validateAll()
        c80, lif = RiskTests.run_calc('breast', pedigree, calcs=['remaining_lifetime', "lifetime"], cwd=self.cwd)
        self.assertEqual(c80, lif)

        # 2. ovarian cancer risk
        pedigree = deepcopy(self.canrisk_pedigree)
        pedigree.validateAll()
        c80, lif = RiskTests.run_calc('ovarian', pedigree, calcs=['remaining_lifetime', "lifetime"], cwd=self.cwd)
        self.assertEqual(c80, lif)

    @pytest.mark.req_WS_VALIDATION_247
    def test_10yr_range_risk(self):
        """ Test ten year range (40-49) risk calculation matches the risk to 50 for a 40 year old. """
        pedigree = deepcopy(self.pedigree)
        pedigree.validateAll()
        calcs1 = Predictions(pedigree, cwd=self.cwd)

        t = pedigree.get_target()
        t.age = '40'
        calcs2 = Predictions(pedigree, cwd=self.cwd)

        c50 = None
        for c in calcs2.cancer_risks:
            if c.get('age') == 50:
                c50 = c

        self.assertEqual(calcs1.ten_yr_cancer_risk[0]['breast cancer risk']['decimal'],
                         c50['breast cancer risk']['decimal'])

    @pytest.mark.req_WS_VALIDATION_248
    def test_incidence_rates(self):
        """ Test prediction of cancer risk and mutation probability for different incidence rates. """
        pedigree = deepcopy(self.pedigree)
        pedigree.validateAll()
        target = pedigree.get_target()
        target.age = 78
        crates = settings.BC_MODEL['CANCER_RATES']

        for cancer_rates in crates.values():
            params = ModelParams(cancer_rates=cancer_rates)
            calcs = Predictions(pedigree, cwd=self.cwd, model_params=params)

            # each gene should have a mutation probability plus a result for no mutations
            self.assertEqual(len(calcs.mutation_probabilties), len(settings.BC_MODEL['GENES']) + 1)

            # risks calculated at different ages:
            self.assertTrue([c.get('age') for c in calcs.cancer_risks] ==
                            [79, 80])

    @pytest.mark.req_WS_VALIDATION_249
    def test_affected_unknown(self):
        """ Test including affected unknown for mother of target to show it increases breast cancer risk. """
        pedigree = deepcopy(self.pedigree)
        target = pedigree.get_target()
        mother = pedigree.get_person(target.mothid)
        mother.yob = str(self.year-55)
        mother.age = "55"
        calcs1 = Predictions(pedigree, cwd=self.cwd)

        def get_c80(calcs):
            for c in calcs.cancer_risks:
                if c.get('age') == 80:
                    return c['breast cancer risk']['decimal']
            return None

        # add affected unknown to mother
        mother.cancers = Cancers(bc1=Cancer("AU"), bc2=Cancer(), oc=Cancer(), prc=Cancer(), pac=Cancer())
        calcs2 = Predictions(pedigree, cwd=self.cwd)
        self.assertGreater(get_c80(calcs2), get_c80(calcs1), 'Mother affected unknown increases BC risk in target')

    @pytest.mark.req_WS_VALIDATION_250
    def test_mutation_frequencies(self):
        """ Test prediction of cancer risk and PV probability for different PV frequencies. """
        pedigree = deepcopy(self.pedigree)
        pedigree.validateAll()
        target = pedigree.get_target()
        target.age = 78
        mutation_frequencies = settings.BC_MODEL['MUTATION_FREQUENCIES']

        for mf in mutation_frequencies.keys():
            if mf == 'Custom':
                continue
            params = ModelParams(population=mf, mutation_frequency=mutation_frequencies[mf])
            calcs = Predictions(pedigree, cwd=self.cwd, model_params=params)

            # each gene should have a mutation probability plus a result for no mutations
            self.assertEqual(len(calcs.mutation_probabilties), len(settings.BC_MODEL['GENES']) + 1)

            # risks calculated at different ages:
            self.assertTrue([c.get('age') for c in calcs.cancer_risks] ==
                            [79, 80])

    @pytest.mark.req_WS_VALIDATION_251
    def test_subproces_err(self):
        """ Test subprocess raises an error when the fortran can not be run. """
        pedigree = deepcopy(self.pedigree)
        pedigree.validateAll()
        with self.assertRaises(FileNotFoundError):
            model_settings = copy.deepcopy(settings.BC_MODEL)
            model_settings['HOME'] = 'xyz'
            Predictions(pedigree, cwd=self.cwd, model_settings=model_settings)

    @pytest.mark.req_WS_VALIDATION_252
    def test_niceness(self):
        """ Test niceness level for pedigree with sibling and large pedigree. """
        pedigree = deepcopy(self.pedigree)
        self.assertEqual(Predictions._get_niceness(pedigree, factor=1), len(pedigree.people))
        self.assertEqual(Predictions._get_niceness(pedigree, factor=0.01), 19)

        # sister
        target = pedigree.get_target()
        sister = Female("FAM1", "F01", "0011", target.fathid, target.mothid, age="22",
                        yob=str(self.year-23), cancers=Cancers(bc1=Cancer("20")))
        pedigree.people.append(sister)
        self.assertEqual(Predictions._get_niceness(pedigree), 1)


class ColumnIdxTests(TestCase):
    ''' Tests for BwaPedigree.get_column_idx and CanRiskPedigree.get_column_idx. '''

    @pytest.mark.req_WS_VALIDATION_190
    def test_bwa_exact_match(self):
        ''' get_column_idx returns the correct index for an exact column name match. '''
        self.assertEqual(BwaPedigree.get_column_idx('Age'), 9)
        self.assertEqual(BwaPedigree.get_column_idx('FamID'), 0)
        self.assertEqual(BwaPedigree.get_column_idx('BRCA1t'), 17)

    @pytest.mark.req_WS_VALIDATION_190
    def test_bwa_case_insensitive(self):
        ''' get_column_idx matches column names case-insensitively. '''
        self.assertEqual(BwaPedigree.get_column_idx('age'), 9)
        self.assertEqual(BwaPedigree.get_column_idx('famid'), 0)
        self.assertEqual(BwaPedigree.get_column_idx('brca1t'), 17)

    @pytest.mark.req_WS_VALIDATION_190
    def test_bwa_not_found(self):
        ''' get_column_idx returns -1 for unknown column names. '''
        self.assertEqual(BwaPedigree.get_column_idx('nonexistent'), -1)

    @pytest.mark.req_WS_VALIDATION_190
    def test_canrisk_exact_match(self):
        ''' CanRiskPedigree.get_column_idx returns the correct index for an exact match. '''
        self.assertEqual(CanRiskPedigree.get_column_idx('Age'), 9)
        self.assertEqual(CanRiskPedigree.get_column_idx('FamID'), 0)

    @pytest.mark.req_WS_VALIDATION_190
    def test_canrisk_case_insensitive(self):
        ''' CanRiskPedigree.get_column_idx matches column names case-insensitively. '''
        self.assertEqual(CanRiskPedigree.get_column_idx('age'), 9)
        self.assertEqual(CanRiskPedigree.get_column_idx('famid'), 0)

    @pytest.mark.req_WS_VALIDATION_190
    def test_canrisk_not_found(self):
        ''' CanRiskPedigree.get_column_idx returns -1 for unknown column names. '''
        self.assertEqual(CanRiskPedigree.get_column_idx('nonexistent'), -1)
