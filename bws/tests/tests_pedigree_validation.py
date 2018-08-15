""" BOADICEA pedigree validation testing.  """
import os
from copy import deepcopy
from datetime import date

from django.test import TestCase

from bws.exceptions import PathologyError, PedigreeError, GeneticTestError,\
    CancerError, PersonError, PedigreeFileError
from bws.pedigree import PedigreeFile, Male, Female
from django.test.utils import override_settings
from bws.cancer import GeneticTest, PathologyTest, PathologyTests, BWSGeneticTests
from bws import cancer
import copy
import re


class ErrorTests(object):
    TEST_BASE_DIR = os.path.dirname(os.path.dirname(__file__))
    TEST_DATA_DIR = os.path.join(TEST_BASE_DIR, 'tests', 'data')

    def setUpErrorTests(self):
        ''' Read in pedigree data. '''
        with open(os.path.join(ErrorTests.TEST_DATA_DIR, "pedigree_data.txt"), "r") as f:
            self.pedigree_data = f.read()
        f.close()
        self.pedigree_file = PedigreeFile(self.pedigree_data)


class PedigreeFileTests(TestCase, ErrorTests):
    """ Tests related to individuals in the pedigree. """

    def setUp(self):
        ''' Read in pedigree data. '''
        super().setUpErrorTests()

    def test_header(self):
        ''' Test pedigree file header. '''
        pedigree_data = copy.copy(self.pedigree_data)
        pedigree_data = pedigree_data.replace('BOADICEA import', 'BOADICEAa import', 1)
        with self.assertRaisesRegex(PedigreeFileError, r"header record in the pedigree file has unexpected characters"):
            PedigreeFile(pedigree_data)

    def test_header2(self):
        ''' Test pedigree file column header. '''
        pedigree_data = copy.copy(self.pedigree_data)
        pedigree_data = pedigree_data.replace('FamID', 'FamIDa', 1)
        with self.assertRaisesRegex(PedigreeFileError, r"headers in the pedigree file contains unexpected characters"):
            PedigreeFile(pedigree_data)

    def test_multiple_pedigrees(self):
        ''' Test multiple pedigrees in a single file. '''
        pedigree_data1 = copy.copy(self.pedigree_data)
        pedigree_data2 = '\n'.join(pedigree_data1.split('\n')[2:]).replace("XXX", "YYY")
        pedigree_data3 = pedigree_data2.replace("YYY", "ZZZ")
        pedigree_data = pedigree_data1 + '\n\n' + pedigree_data2 + '\n\n' + pedigree_data3
        pedigrees = PedigreeFile(pedigree_data).pedigrees
        self.assertEqual(len(pedigrees), 3)
        for p in pedigrees:
            warnings = PedigreeFile.validate(p)
            self.assertEqual(len(warnings), 0)

    def test_columns(self):
        ''' Test multiple pedigrees in a single file. '''
        pedigree_data = copy.copy(self.pedigree_data)
        pedigree_data = pedigree_data.replace('F1', 'F1    F1')
        with self.assertRaisesRegex(PedigreeFileError, r"data record has an unexpected number of data items"):
            PedigreeFile(pedigree_data)


class PersonTests(TestCase, ErrorTests):
    """ Tests related to individuals in the pedigree. """

    def setUp(self):
        ''' Read in pedigree data. '''
        super().setUpErrorTests()

    def test_name(self):
        """ Test an error is raised if the individuals name is not an alphanumeric string. """
        pedigree_file = deepcopy(self.pedigree_file)
        m2 = pedigree_file.pedigrees[0].get_person_by_name('M2')
        m2.name = "M<2"
        with self.assertRaisesRegex(PersonError, r"is unspecified or is not an alphanumeric string"):
            PedigreeFile.validate(pedigree_file.pedigrees)

    def test_indvid_id(self):
        """ Test an error is raised if the individuals ID is not an alphanumeric or is too large.  """
        pedigree_file = deepcopy(self.pedigree_file)
        pedigree = pedigree_file.pedigrees[0]
        m2 = pedigree.get_person_by_name('M2')
        m2.pid = "XXXXXXXX"
        with self.assertRaisesRegex(PersonError, r"Individual identifiers must be alphanumeric strings with a max"):
            m2.validate(pedigree)

        m2.pid = ""
        with self.assertRaisesRegex(PersonError, r"Individual identifiers must be alphanumeric strings with a max"):
            m2.validate(pedigree)

        m2.pid = "000"
        with self.assertRaisesRegex(PersonError, r"Individual identifiers must be alphanumeric strings with a max"):
            m2.validate(pedigree)

        m2.pid = "X*A"
        with self.assertRaisesRegex(PersonError, r"Individual identifiers must be alphanumeric strings with a max"):
            m2.validate(pedigree)

    def test_fathid_id(self):
        """ Test an error is raised if the father ID is not an alphanumeric or is too large.  """
        pedigree_file = deepcopy(self.pedigree_file)
        pedigree = pedigree_file.pedigrees[0]
        m2 = pedigree.get_person_by_name('M2')
        m2.fathid = "XXXXXXXX"
        with self.assertRaisesRegex(PersonError, r"Father identifier .* unexpected character"):
            m2.validate(pedigree)

        m2.fathid = ""
        with self.assertRaisesRegex(PersonError, r"Father identifier .* unexpected character"):
            m2.validate(pedigree)

        m2.fathid = "X*A"
        with self.assertRaisesRegex(PersonError, r"Father identifier .* unexpected character"):
            m2.validate(pedigree)

    def test_mothid_id(self):
        """ Test an error is raised if the mother ID is not an alphanumeric or is too large.  """
        pedigree_file = deepcopy(self.pedigree_file)
        pedigree = pedigree_file.pedigrees[0]
        m2 = pedigree.get_person_by_name('M2')
        m2.mothid = "XXXXXXXX"
        with self.assertRaisesRegex(PersonError, r"Mother identifier .* unexpected character"):
            m2.validate(pedigree)

        m2.mothid = ""
        with self.assertRaisesRegex(PersonError, r"Mother identifier .* unexpected character"):
            m2.validate(pedigree)

        m2.mothid = "X*A"
        with self.assertRaisesRegex(PersonError, r"Mother identifier .* unexpected character"):
            m2.validate(pedigree)

    def test_parent_unspecified(self):
        """ Test an error is raised if only one of the parents is specified. """
        pedigree_file = deepcopy(self.pedigree_file)
        pedigree = pedigree_file.pedigrees[0]
        m2 = pedigree.get_person_by_name('F1')
        m2.mothid = "0"
        with self.assertRaisesRegex(PersonError, r"only one parent specified"):
            m2.validate(pedigree)

    def test_father_sex(self):
        ''' Test an error is raised if the sex of the father is not 'M'. '''
        # change sex of father to 'F'
        pedigree_file = deepcopy(self.pedigree_file)
        m2 = pedigree_file.pedigrees[0].get_person_by_name('M2')
        pedigree_file.pedigrees[0].people.remove(m2)
        pedigree_file.pedigrees[0].people.append(Female(m2.famid, m2.name, m2.pid, m2.fathid, m2.mothid))
        with self.assertRaisesRegex(PersonError, r"All fathers in the pedigree must have sex specified as 'M'"):
            PedigreeFile.validate(pedigree_file.pedigrees)

    def test_mother_sex(self):
        ''' Test an error is raised if the sex of the mother is not 'F'. '''
        # change sex of mother to 'M'
        pedigree_file = deepcopy(self.pedigree_file)
        f2 = pedigree_file.pedigrees[0].get_person_by_name('F2')
        pedigree_file.pedigrees[0].people.remove(f2)
        pedigree_file.pedigrees[0].people.append(Male(f2.famid, f2.name, f2.pid, f2.fathid, f2.mothid))
        with self.assertRaisesRegex(PersonError, r"All mothers in the pedigree must have sex specified as 'F'"):
            PedigreeFile.validate(pedigree_file.pedigrees)

    def test_missing_mother(self):
        ''' Test an error is raised if the mother is missing. '''
        # remove mother F2 from pedigree
        pedigree_file = deepcopy(self.pedigree_file)
        f2 = pedigree_file.pedigrees[0].get_person_by_name('F2')
        pedigree_file.pedigrees[0].people.remove(f2)
        with self.assertRaisesRegex(PersonError, r"The mother (.*) is missing"):
            PedigreeFile.validate(pedigree_file.pedigrees)

    def test_missing_father(self):
        ''' Test an error is raised if the father is missing. '''
        # remove father M2 from pedigree
        pedigree_file = deepcopy(self.pedigree_file)
        PedigreeFile.validate(pedigree_file.pedigrees)
        m2 = pedigree_file.pedigrees[0].get_person_by_name('M2')
        pedigree_file.pedigrees[0].people.remove(m2)
        with self.assertRaisesRegex(PersonError, r"The father (.*) is missing"):
            PedigreeFile.validate(pedigree_file.pedigrees)

    def test_dead(self):
        ''' Test an error is raised if the dead attribute is incorrectly specified (not 0 or 1). '''
        pedigree_file = deepcopy(self.pedigree_file)
        PedigreeFile.validate(pedigree_file.pedigrees)
        m2 = pedigree_file.pedigrees[0].get_person_by_name('M2')
        m2.dead = "3"
        with self.assertRaisesRegex(PersonError, r"alive must be specified as '0', and dead specified as '1'"):
            PedigreeFile.validate(pedigree_file.pedigrees)

    def test_age_yob(self):
        ''' Test a warning is returned if the yob and age attribute is not specified. '''
        pedigree_file = deepcopy(self.pedigree_file)
        m2 = pedigree_file.pedigrees[0].get_person_by_name('M2')
        m2.age = "0"
        m2.cancers.diagnoses.prc.age = "0"
        warnings = PedigreeFile.validate(pedigree_file.pedigrees)
        self.assertRegex(warnings[0], "year of birth and age at last follow up must be specified")

    def test_age(self):
        ''' Test an error is raised if the age attribute is incorrectly specified (not 0 or 1-MAX_AGE). '''
        pedigree_file = deepcopy(self.pedigree_file)
        PedigreeFile.validate(pedigree_file.pedigrees)
        m2 = pedigree_file.pedigrees[0].get_person_by_name('M2')
        m2.age = "222"
        with self.assertRaisesRegex(PersonError, r"Ages must be specified with as '0' for unknown, or in the range"):
            PedigreeFile.validate(pedigree_file.pedigrees)

    def test_yob(self):
        ''' Test error is raised if the year of birth attribute is incorrectly specified. '''
        pedigree_file = deepcopy(self.pedigree_file)
        PedigreeFile.validate(pedigree_file.pedigrees)
        m2 = pedigree_file.pedigrees[0].get_person_by_name('M2')
        m2.yob = "1200"
        with self.assertRaisesRegex(PersonError, r"Years of birth must be in the range"):
            PedigreeFile.validate(pedigree_file.pedigrees)

        m2.yob = "xyz"
        with self.assertRaisesRegex(PersonError, r"Years of birth must be in the range"):
            PedigreeFile.validate(pedigree_file.pedigrees)

        m2.yob = str(date.today().year+1)
        with self.assertRaisesRegex(PersonError, r"Years of birth must be in the range"):
            PedigreeFile.validate(pedigree_file.pedigrees)

    def test_ashkn(self):
        """ Test an error is raised if the Ashkenazi origin is incorrectly assigned (i.e. not 0 or 1). """
        pedigree_file = deepcopy(self.pedigree_file)
        f1 = pedigree_file.pedigrees[0].get_person_by_name('F1')
        f1.ashkn = "x"
        with self.assertRaisesRegex(PersonError, r"invalid Ashkenazi"):
            PedigreeFile.validate(pedigree_file.pedigrees)

    @override_settings(MAX_NUMBER_OF_SIBS_PER_NUCLEAR_FAMILY=1)
    def test_siblings(self):
        """ Check max no. of siblings not exceeded. """
        pedigree_file = deepcopy(self.pedigree_file)
        pedigree = pedigree_file.pedigrees[0]
        f1 = pedigree.get_person_by_name('F1')
        self.assertEqual(len(pedigree.get_siblings(f1)[0]), 0)
        pedigree.people.append(Male(f1.famid, "M1A", "111", f1.fathid, f1.mothid))
        pedigree.people.append(Male(f1.famid, "M1B", "112", f1.fathid, f1.mothid))
        self.assertEqual(len(pedigree.get_siblings(f1)[0]), 2)
        with self.assertRaisesRegex(PersonError, r"exeeded the maximum number of siblings"):
            PedigreeFile.validate(pedigree_file.pedigrees)

    @override_settings(MAX_NUMBER_OF_SIBS_PER_NUCLEAR_FAMILY_WITH_SAME_YOB=1)
    def test_siblings_same_yob(self):
        """ Check max no. of siblings with same year of birth is not exceeded. """
        pedigree_file = deepcopy(self.pedigree_file)
        pedigree = pedigree_file.pedigrees[0]
        f1 = pedigree.get_person_by_name('F1')
        self.assertEqual(len(pedigree.get_siblings(f1)[1]), 0)
        pedigree.people.append(Male(f1.famid, "M1A", "111", f1.fathid, f1.mothid, yob=f1.yob))
        pedigree.people.append(Male(f1.famid, "M1B", "112", f1.fathid, f1.mothid, yob=f1.yob))
        self.assertEqual(len(pedigree.get_siblings(f1)[1]), 2)
        with self.assertRaisesRegex(PersonError, r"siblings with the same year of birth exceeded"):
            PedigreeFile.validate(pedigree_file.pedigrees)


class PedigreeTests(TestCase, ErrorTests):
    """ Tests related to pedigree input data. """

    def setUp(self):
        ''' Read in pedigree data. '''
        super().setUpErrorTests()

    def test_target(self):
        """ Test an error is raised if target column is not between 0 and 1. """
        pedigree_data = copy.copy(self.pedigree_data)
        pedigree_data = re.sub(r"F1\s+1", 'F1    3', pedigree_data)
        with self.assertRaisesRegex(PedigreeError, r"value in the Target data column"):
            PedigreeFile(pedigree_data)

    def test_no_target(self):
        """ Test an error is raised if there is no target. """
        pedigree_data = copy.copy(self.pedigree_data)
        pedigree_data = re.sub(r"F1\s+1", 'F1    0', pedigree_data)
        with self.assertRaisesRegex(PedigreeError, r"has either no index or more than 1"):
            PedigreeFile(pedigree_data)

    def test_duplicate_person(self):
        """ Test an error is raised if an individual appears more than once. """
        pedigree_data1 = copy.copy(self.pedigree_data)
        pedigree_data2 = '\n'.join(pedigree_data1.split('\n')[2:3])
        pedigree_data1 = pedigree_data1 + pedigree_data2
        with self.assertRaisesRegex(PedigreeError, r"appears more than once in the pedigree"):
            PedigreeFile(pedigree_data1)

    @override_settings(MAX_PEDIGREE_SIZE=2)
    def test_no_family_members(self):
        """ Test an error is raised if number of people in a pedigree is > MAX_PEDIGREE_SIZE. """
        with self.assertRaisesRegex(PedigreeError, r"unexpected number of family members"):
            PedigreeFile(self.pedigree_data)

    def test_sex(self):
        """ Test an error is raised if individuals sex is not M or F. """
        pedigree_data = copy.copy(self.pedigree_data)
        pedigree_data = re.sub(r"\s+F\s+", '   FF   ', pedigree_data)
        with self.assertRaisesRegex(PedigreeError, r"individuals sex must be specified as"):
            PedigreeFile(pedigree_data)

    def test_famid(self):
        """ Test an error is raised if the family ID is incorrectly specified. """
        pedigree_file = deepcopy(self.pedigree_file)
        pedigree_file.pedigrees[0].famid = "XXXXXXXXXXXXXX"
        with self.assertRaisesRegex(PedigreeError, r"Family IDs must be specified with between 1 and "):
            PedigreeFile.validate(pedigree_file.pedigrees)

        pedigree_file.pedigrees[0].famid = "X?"
        with self.assertRaisesRegex(PedigreeError, r"Family IDs must be specified with between 1 and "):
            PedigreeFile.validate(pedigree_file.pedigrees)

        pedigree_file.pedigrees[0].famid = "00"
        with self.assertRaisesRegex(PedigreeError, r"Family IDs must be specified with between 1 and "):
            PedigreeFile.validate(pedigree_file.pedigrees)

    def test_target_yob(self):
        """ The target must be assigned a valid year of birth. """
        pedigree_file = deepcopy(self.pedigree_file)
        f1 = pedigree_file.pedigrees[0].get_person_by_name('F1')
        f1.yob = "0"
        with self.assertRaisesRegex(PedigreeError, r"This person must be assigned a valid year of birth"):
            PedigreeFile.validate(pedigree_file.pedigrees)

    def test_target_age(self):
        """ The target must be assigned a valid year of birth. """
        pedigree_file = deepcopy(self.pedigree_file)
        f1 = pedigree_file.pedigrees[0].get_person_by_name('F1')
        f1.age = "0"
        with self.assertRaisesRegex(PedigreeError, r"This person must be assigned an age"):
            PedigreeFile.validate(pedigree_file.pedigrees)

    def test_target_risk_prob_calcs(self):
        """ Check if an error is raised if the carrier probabilities and cancer risk can not be
        calculated.

        If the target has a positive genetic test carrier probs cannot be calculated. Had bilateral
        BC or OC or PanC, or she is too old. """
        pedigree_file = deepcopy(self.pedigree_file)
        f1 = pedigree_file.pedigrees[0].get_person_by_name('F1')
        f1.gtests.brca1.result = "P"
        f1.cancers.diagnoses.bc2.age = 18
        with self.assertRaisesRegex(PedigreeError, r"cannot compute mutation carrier probabilities"):
            PedigreeFile.validate(pedigree_file.pedigrees)

    def test_ununconnected(self):
        pedigree_file = deepcopy(self.pedigree_file)
        m2 = pedigree_file.pedigrees[0].get_person_by_name('M2')
        pedigree_file.pedigrees[0].people.append(Male(m2.famid, "M1A", "111", "0", "0"))
        with self.assertRaisesRegex(PedigreeError, r"family members are not physically connected to the target"):
            PedigreeFile.validate(pedigree_file.pedigrees)

    def test_mztwin_pair(self):
        """ Check if an error is raised if the number of people specified in a set of twins is not 2. """
        pedigree_file = deepcopy(self.pedigree_file)
        m2 = pedigree_file.pedigrees[0].get_person_by_name('M2')
        m2.mztwin = "1"
        with self.assertRaisesRegex(PedigreeError, r"Only MZ twins are permitted in the pedigree"):
            PedigreeFile.validate(pedigree_file.pedigrees)

    @override_settings(MAX_NUMBER_MZ_TWIN_PAIRS=0)
    def test_mztwin_number(self):
        """ Check if an error is raised if the number of people specified in a set of twins is not 2. """
        pedigree_file = deepcopy(self.pedigree_file)
        self._add_twin(pedigree_file.pedigrees[0])
        with self.assertRaisesRegex(PedigreeError, r"Maximum number of MZ twin pairs has been exceeded"):
            PedigreeFile.validate(pedigree_file.pedigrees)

    def test_mztwin_key(self):
        """ Check if an error is raised if MZ twin characters are not valid. """
        pedigree_file = deepcopy(self.pedigree_file)
        self._add_twin(pedigree_file.pedigrees[0], twin_key="X")
        with self.assertRaisesRegex(PedigreeError, r"MZ twins must be identified using one"):
            PedigreeFile.validate(pedigree_file.pedigrees)

    def test_mztwin_parents(self):
        """ Check if an error is raised if MZ twin have different parents. """
        pedigree_file = deepcopy(self.pedigree_file)
        pedigree = pedigree_file.pedigrees[0]
        m2 = pedigree.get_person_by_name('M2')
        self._add_twin(pedigree, twin=m2)
        m2.mothid = "313"
        pedigree.people.append(Female(m2.famid, "F3A", "313", "0", "0"))
        with self.assertRaisesRegex(PedigreeError, r"MZ twins must have the same parents"):
            PedigreeFile.validate(pedigree_file.pedigrees)

    def test_mztwin_yob(self):
        """ Check if an error is raised if MZ twin have different years of birth. """
        pedigree_file = deepcopy(self.pedigree_file)
        pedigree = pedigree_file.pedigrees[0]
        m2 = pedigree.get_person_by_name('M2')
        self._add_twin(pedigree, twin=m2)
        m2.yob = str(int(m2.yob) + 2)
        with self.assertRaisesRegex(PedigreeError, r"MZ twins must have the same year of birth"):
            PedigreeFile.validate(pedigree_file.pedigrees)

    def test_mztwin_age(self):
        """ Check if an error is raised if MZ twin are alive and have different ages at last followup. """
        pedigree_file = deepcopy(self.pedigree_file)
        pedigree = pedigree_file.pedigrees[0]
        m2 = pedigree.get_person_by_name('M2')
        self._add_twin(pedigree, twin=m2)
        m2.age = "44"
        with self.assertRaisesRegex(PedigreeError, r"If both MZ twins are alive, they must have the same age at last"):
            PedigreeFile.validate(pedigree_file.pedigrees)

    def test_mztwin_sex(self):
        """ Check if an error is raised if MZ twins are not the same sex. """
        pedigree_file = deepcopy(self.pedigree_file)
        pedigree = pedigree_file.pedigrees[0]
        m2 = pedigree.get_person_by_name('M2')
        self._add_twin(pedigree, twin=m2)

        pedigree_file.pedigrees[0].people.remove(m2)
        pedigree_file.pedigrees[0].people.append(Female(m2.famid, m2.name, m2.pid, m2.fathid, m2.mothid,
                                                        mztwin=m2.mztwin, age=m2.age, yob=m2.yob,
                                                        cancers=m2.cancers))
        with self.assertRaisesRegex(PedigreeError, r"MZ twins must have the same sex"):
            PedigreeFile.validate(pedigree_file.pedigrees)

    def test_mztwin_genetic_tests(self):
        """ Check if an error is raised if MZ twins have different genetic test results. """
        pedigree_file = deepcopy(self.pedigree_file)
        pedigree = pedigree_file.pedigrees[0]
        m2 = pedigree.get_person_by_name('M2')

        m2.gtests = BWSGeneticTests._make([GeneticTest("S", "P")
                                           if i == 0 else GeneticTest("0", "0")
                                           for i in range(len(BWSGeneticTests._fields))])
        self._add_twin(pedigree, twin=m2)
        m2.gtests = BWSGeneticTests._make([GeneticTest("S", "N")
                                           if i == 0 else GeneticTest("0", "0")
                                           for i in range(len(BWSGeneticTests._fields))])

        with self.assertRaisesRegex(PedigreeError, r"genetic test results for these individuals are different"):
            PedigreeFile.validate(pedigree_file.pedigrees)

    def _add_twin(self, pedigree, twin=None, twin_key="1"):
        """ Utility for creating a twin given a person in a pedigree.  """
        if twin is None:
            twin = pedigree.get_person_by_name('M2')
        twin.mztwin = twin_key
        twin.fathid = "311"
        twin.mothid = "312"
        pedigree.people.append(Male(twin.famid, "M3", "311", "0", "0"))
        pedigree.people.append(Female(twin.famid, "F3", "312", "0", "0"))
        pedigree.people.append(Male(twin.famid, "M2A", "111", "311", "312", mztwin=twin_key,
                                    age=twin.age, yob=twin.yob, gtests=deepcopy(twin.gtests),
                                    cancers=deepcopy(twin.cancers)))


class CancerTests(TestCase, ErrorTests):

    def setUp(self):
        ''' Read in pedigree data. '''
        super().setUpErrorTests()

    def test_cancer_age(self):
        ''' Test an error is raised if the cancer diagnosis age is incorrectly specified (not 0 or 1-MAX_AGE). '''
        pedigree_file = deepcopy(self.pedigree_file)
        f1 = pedigree_file.pedigrees[0].get_person_by_name('F1')
        f1.cancers.diagnoses.bc1.age = "xyz"
        with self.assertRaisesRegex(CancerError, r"Age at cancer diagnosis"):
            PedigreeFile.validate(pedigree_file.pedigrees)

        f1.cancers.diagnoses.bc1.age = "199"
        with self.assertRaisesRegex(CancerError, r"Age at cancer diagnosis"):
            PedigreeFile.validate(pedigree_file.pedigrees)

    def test_cancer_age2(self):
        ''' Test an error is raised if the cancer diagnosis age is greater than age at last follow up. '''
        pedigree_file = deepcopy(self.pedigree_file)
        m2 = pedigree_file.pedigrees[0].get_person_by_name('M2')
        m2.cancers.diagnoses.bc1.age = "59"
        m2.age = "53"
        with self.assertRaisesRegex(CancerError, r"diagnosis that exceeds age at last follow up"):
            PedigreeFile.validate(pedigree_file.pedigrees)

    def test_cancer_yob(self):
        ''' Test an error is raised if the yob is not specified when a cancer is diagnosed. '''
        pedigree_file = deepcopy(self.pedigree_file)
        m2 = pedigree_file.pedigrees[0].get_person_by_name('M2')
        m2.yob = "0"
        with self.assertRaisesRegex(CancerError, r"diagnosed with cancer but has no year of birth specified"):
            PedigreeFile.validate(pedigree_file.pedigrees)

    def test_male_with_ovarian_cancer(self):
        """ Test an error is raised if a male has been designated ovarian cancer. """
        pedigree_file = deepcopy(self.pedigree_file)
        m2 = pedigree_file.pedigrees[0].get_person_by_name('M2')
        m2.cancers.diagnoses.oc.age = "45"
        with self.assertRaisesRegex(CancerError, r"male but has been assigned an ovarian"):
            PedigreeFile.validate(pedigree_file.pedigrees)

    def test_female_with_prostate_cancer(self):
        """ Test an error is raised if a female has been designated prostate cancer. """
        pedigree_file = deepcopy(self.pedigree_file)
        f1 = pedigree_file.pedigrees[0].get_person_by_name('F1')
        f1.cancers.diagnoses.prc.age = "21"
        with self.assertRaisesRegex(CancerError, r"female but has been assigned an prostate"):
            PedigreeFile.validate(pedigree_file.pedigrees)

    def test_bc1_missing(self):
        """ Test an error is raised if the age of a second breast cancer is present but age of first is missing. """
        pedigree_file = deepcopy(self.pedigree_file)
        f1 = pedigree_file.pedigrees[0].get_person_by_name('F1')
        f1.cancers.diagnoses.bc1.age = "0"
        f1.cancers.diagnoses.bc2.age = "21"
        with self.assertRaisesRegex(CancerError, r"contralateral breast cancer, (.*) " +
                                    "first breast cancer is missing"):
            PedigreeFile.validate(pedigree_file.pedigrees)

    def test_bc1_missing2(self):
        """ Test an error is raised if the age of a second breast cancer is AU but age of first is missing. """
        pedigree_file = deepcopy(self.pedigree_file)
        f1 = pedigree_file.pedigrees[0].get_person_by_name('F1')
        f1.cancers.diagnoses.bc1.age = "0"
        f1.cancers.diagnoses.bc2.age = "AU"
        with self.assertRaisesRegex(CancerError, r"contralateral breast cancer, (.*) " +
                                    "first breast cancer is missing"):
            PedigreeFile.validate(pedigree_file.pedigrees)

    def test_bc1_and_bc2(self):
        """ Test an error is raised if the age of a first breast cancer exceeds that of the second. """
        pedigree_file = deepcopy(self.pedigree_file)
        f1 = pedigree_file.pedigrees[0].get_person_by_name('F1')
        f1.cancers.diagnoses.bc1.age = "22"
        f1.cancers.diagnoses.bc2.age = "21"
        with self.assertRaisesRegex(CancerError, r"contralateral breast cancer, (.*) " +
                                    "age at diagnosis of the first breast cancer exceeds"):
            PedigreeFile.validate(pedigree_file.pedigrees)


class GeneticTestTests(TestCase, ErrorTests):
    """ Test related to genetic test results. """

    def setUp(self):
        ''' Read in pedigree data. '''
        super().setUpErrorTests()

    def test_type(self):
        """ Check that the genetic test type is valid. """
        pedigree_file = deepcopy(self.pedigree_file)
        f1 = pedigree_file.pedigrees[0].get_person_by_name('F1')
        f1.gtests.brca1.test_type = "X"
        with self.assertRaisesRegex(GeneticTestError, "invalid genetic test type"):
            PedigreeFile.validate(pedigree_file.pedigrees)

    def test_type_specified(self):
        """ Check if there is a genetic test result check the test type is specified. """
        pedigree_file = deepcopy(self.pedigree_file)
        f1 = pedigree_file.pedigrees[0].get_person_by_name('F1')
        f1.gtests.brca1.test_type = "0"
        f1.gtests.brca1.result = "P"
        with self.assertRaisesRegex(GeneticTestError, "genetic test type has not been specified"):
            PedigreeFile.validate(pedigree_file.pedigrees)

    def test_result(self):
        """ Check that the mutation status is valid. """
        pedigree_file = deepcopy(self.pedigree_file)
        f1 = pedigree_file.pedigrees[0].get_person_by_name('F1')
        f1.gtests.brca1.test_type = "S"
        f1.gtests.brca1.result = "X"
        with self.assertRaisesRegex(GeneticTestError, "invalid genetic test result"):
            PedigreeFile.validate(pedigree_file.pedigrees)

    def test_result_specified(self):
        """ Check that the mutation status is specified if tested. """
        pedigree_file = deepcopy(self.pedigree_file)
        f1 = pedigree_file.pedigrees[0].get_person_by_name('F1')
        f1.gtests.brca1.test_type = "S"
        f1.gtests.brca1.result = "0"
        with self.assertRaisesRegex(GeneticTestError, "corresponding test result has not been specified"):
            PedigreeFile.validate(pedigree_file.pedigrees)


class PathologyTestTests(TestCase, ErrorTests):
    """ Tests related to pathology test results. """

    def setUp(self):
        ''' Read in pedigree data. '''
        super().setUpErrorTests()

    def test_pathology_status(self):
        """ Check that the pathology results are correctly set (0, N, P). """
        pedigree_file = deepcopy(self.pedigree_file)
        f1 = pedigree_file.pedigrees[0].get_person_by_name('F1')
        f1.pathology.er.result = 'X'
        with self.assertRaisesRegex(PathologyError, "must be 'N' for negative, 'P' for positive, or '0' for unknown"):
            PedigreeFile.validate(pedigree_file.pedigrees)

    def test_patholgy_bc1_setting(self):
        """ Check that pathology test results are only provided for family members with a first breast cancer """
        pedigree_file = deepcopy(self.pedigree_file)
        f1 = pedigree_file.pedigrees[0].get_person_by_name('F1')
        f1.cancers.diagnoses.bc1.age = "0"
        f1.pathology.er.result = "P"
        with self.assertRaisesRegex(PathologyError, "Pathology test results can only be assigned to family " +
                                    "members who have developed breast cancer."):
            PedigreeFile.validate(pedigree_file.pedigrees)

    def test_er_unspecified(self):
        """
        If an individual has breast cancer and ER is unspecified but another pathology parameter has
        been specified, warn that no pathology data will be used.
        """
        pedigree_file = deepcopy(self.pedigree_file)
        f1 = pedigree_file.pedigrees[0].get_person_by_name('F1')
        f1. pathology = PathologyTests(
                    er=PathologyTest(cancer.ESTROGEN_RECEPTOR_TEST, result="0"),
                    pr=PathologyTest(cancer.PROGESTROGEN_RECEPTOR_TEST, result="P"),
                    her2=PathologyTest(cancer.HER2_TEST, result="0"),
                    ck14=PathologyTest(cancer.CK14_TEST, result="0"),
                    ck56=PathologyTest(cancer.CK56_TEST, result="0"))
        warnings = PedigreeFile.validate(pedigree_file.pedigrees)
        self.assertRegex(warnings[0], "this individual's pathology information will not be taken into account")

    def test_er_positive(self):
        """
        If an individual has breast cancer and ER positive, and PR or HER2 or CK14 or CK5/6 specified
        generate a warning
        """
        pedigree_file = deepcopy(self.pedigree_file)
        f1 = pedigree_file.pedigrees[0].get_person_by_name('F1')
        f1. pathology = PathologyTests(
                    er=PathologyTest(cancer.ESTROGEN_RECEPTOR_TEST, result="P"),
                    pr=PathologyTest(cancer.PROGESTROGEN_RECEPTOR_TEST, result="P"),
                    her2=PathologyTest(cancer.HER2_TEST, result="0"),
                    ck14=PathologyTest(cancer.CK14_TEST, result="0"),
                    ck56=PathologyTest(cancer.CK56_TEST, result="0"))
        warnings = PedigreeFile.validate(pedigree_file.pedigrees)
        self.assertRegex(warnings[0], "ER positive, where an additional pathology parameter ")

    def test_pr_or_her_unspecified(self):
        """
        If an individual has breast cancer and ER negative and PR status is specified but HER2 status
        is unspecified (or vice versa) report a warning.
        """
        pedigree_file = deepcopy(self.pedigree_file)
        f1 = pedigree_file.pedigrees[0].get_person_by_name('F1')
        f1. pathology = PathologyTests(
                    er=PathologyTest(cancer.ESTROGEN_RECEPTOR_TEST, result="N"),
                    pr=PathologyTest(cancer.PROGESTROGEN_RECEPTOR_TEST, result="P"),
                    her2=PathologyTest(cancer.HER2_TEST, result="0"),
                    ck14=PathologyTest(cancer.CK14_TEST, result="0"),
                    ck56=PathologyTest(cancer.CK56_TEST, result="0"))
        warnings = PedigreeFile.validate(pedigree_file.pedigrees)
        self.assertRegex(warnings[0], "PR status is specified but HER2 status is unspecified")

    def test_ck14_ck56_unspecified(self):
        """
        If an individual has breast cancer and either CK14 or CK5/6 has been specified (one without the
        other) generate a warning.
        """
        pedigree_file = deepcopy(self.pedigree_file)
        f1 = pedigree_file.pedigrees[0].get_person_by_name('F1')
        f1. pathology = PathologyTests(
                    er=PathologyTest(cancer.ESTROGEN_RECEPTOR_TEST, result="N"),
                    pr=PathologyTest(cancer.PROGESTROGEN_RECEPTOR_TEST, result="N"),
                    her2=PathologyTest(cancer.HER2_TEST, result="N"),
                    ck14=PathologyTest(cancer.CK14_TEST, result="N"),
                    ck56=PathologyTest(cancer.CK56_TEST, result="0"))
        warnings = PedigreeFile.validate(pedigree_file.pedigrees)
        self.assertRegex(warnings[0], "only CK14 or CK5/6 status has been specified")

    def test_triple_negative(self):
        """
        If an individual has breast cancer and not triple Negative (ER, PR, HER2) but
        CK14 and CK5/6 are specified generate a warning.
        """
        pedigree_file = deepcopy(self.pedigree_file)
        f1 = pedigree_file.pedigrees[0].get_person_by_name('F1')
        f1. pathology = PathologyTests(
                    er=PathologyTest(cancer.ESTROGEN_RECEPTOR_TEST, result="N"),
                    pr=PathologyTest(cancer.PROGESTROGEN_RECEPTOR_TEST, result="0"),
                    her2=PathologyTest(cancer.HER2_TEST, result="0"),
                    ck14=PathologyTest(cancer.CK14_TEST, result="N"),
                    ck56=PathologyTest(cancer.CK56_TEST, result="N"))
        warnings = PedigreeFile.validate(pedigree_file.pedigrees)
        self.assertRegex(warnings[0], "CK14 or CK5/6 status is specified but the breast cancer pathology is not triple")
