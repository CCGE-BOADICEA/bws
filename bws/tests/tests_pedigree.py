""" Mutation risk and probability calculation testing. """
from datetime import date
from django.test import TestCase
from bws.pedigree import Female, PedigreeFile, BwaPedigree, CanRiskPedigree
from copy import deepcopy
from bws.calcs import Predictions, ModelParams
from bws.cancer import Cancer, Cancers, CanRiskGeneticTests
from django.conf import settings
import tempfile
import shutil
import os
import filecmp
import copy


class WritePedigree(TestCase):
    TEST_BASE_DIR = os.path.dirname(os.path.dirname(__file__))
    TEST_DATA_DIR = os.path.join(TEST_BASE_DIR, 'tests', 'data')

    def test_write(self):
        """ Write out BOADICEA pedigree file and compare with the original file. """
        with open(os.path.join(WritePedigree.TEST_DATA_DIR, "d3.bwa"), "r") as f1:
            pedigree_data = f1.read()

        pedigree_file = PedigreeFile(pedigree_data)
        pedigree = pedigree_file.pedigrees[0]

        f2 = tempfile.NamedTemporaryFile(mode='w')
        pedigree.write_boadicea_file_header(bwa_file=f2)
        pedigree.write_boadicea_file(bwa_file=f2)
        self.assertTrue(filecmp.cmp(f1.name, f2.name, shallow=False))


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

    def test_calculations(self):
        """ Test prediction of cancer risk and mutation probability. """
        pedigree = deepcopy(self.pedigree)
        PedigreeFile.validate(pedigree)
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

    def test_calculations2(self):
        """ Test prediction of cancer risk and mutation probability. """
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

        PedigreeFile.validate(pedigree)
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

    def test_ovarian_calculations(self):
        """ Test prediction of cancer risk and mutation probability. """
        target = Female("FAM1", "F0", "001", "002", "003", target="1", age="20",
                        yob=str(self.year-20), cancers=Cancers(),
                        gtests=CanRiskGeneticTests.default_factory())
        pedigree = CanRiskPedigree(people=[target])
        # parents
        (_father, _mother) = pedigree.add_parents(target, gtests=CanRiskGeneticTests.default_factory())
        PedigreeFile.validate(pedigree)
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

    def test_lifetime_risk(self):
        """ Test lifetime risk calculation matches the risk to 80 for a 20 year old. """
        # 1. breast cancer risk
        pedigree = deepcopy(self.pedigree)
        t = pedigree.get_target()
        t.age = '20'
        PedigreeFile.validate(pedigree)
        c80, lif = RiskTests.run_calc('breast', pedigree, calcs=['remaining_lifetime', "lifetime"], cwd=self.cwd)
        self.assertEqual(c80, lif)

        # 2. ovarian cancer risk
        pedigree = deepcopy(self.canrisk_pedigree)
        PedigreeFile.validate(pedigree)
        c80, lif = RiskTests.run_calc('ovarian', pedigree, calcs=['remaining_lifetime', "lifetime"], cwd=self.cwd)
        self.assertEqual(c80, lif)

    def test_10yr_range_risk(self):
        """ Test ten year range (40-49) risk calculation matches the risk to 50 for a 40 year old. """
        pedigree = deepcopy(self.pedigree)
        PedigreeFile.validate(pedigree)
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

    def test_incidence_rates(self):
        """ Test prediction of cancer risk and mutation probability for different incidence rates. """
        pedigree = deepcopy(self.pedigree)
        PedigreeFile.validate(pedigree)
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

    def test_mutation_frequencies(self):
        """ Test prediction of cancer risk and mutation probability for different mutation frequencies. """
        pedigree = deepcopy(self.pedigree)
        PedigreeFile.validate(pedigree)
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

    def test_subproces_err(self):
        """ Test subprocess raises an error when the fortran can not be run. """
        pedigree = deepcopy(self.pedigree)
        PedigreeFile.validate(pedigree)
        with self.assertRaises(FileNotFoundError):
            model_settings = copy.deepcopy(settings.BC_MODEL)
            model_settings['HOME'] = 'xyz'
            Predictions(pedigree, cwd=self.cwd, model_settings=model_settings)

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
