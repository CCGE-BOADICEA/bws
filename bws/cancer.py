""" Cancer, pathology and genetic testing """
import collections
import re
from bws.exceptions import GeneticTestError, PedigreeFileError,\
    PathologyError, CancerError
from django.conf import settings


REGEX_PATHOLOGY_TEST_OPTION = re.compile("^([1-5])$")
REGEX_PATHOLOGY_STATUS = re.compile("^[0NP]$")

REGEX_GENETIC_TEST_TYPE = re.compile("^[0ST]$")
REGEX_BOADICEA_FORMAT_4_GENETIC_TEST_RESULT = re.compile("^[0NP]$")
REGEX_GENETIC_TEST_TYPE_IS_TESTED = re.compile("^[ST]$")

# pathology tests
ESTROGEN_RECEPTOR_TEST = "1"      # Estrogen Receptor test
PROGESTROGEN_RECEPTOR_TEST = "2"  # Progestrogen Receptor test
HER2_TEST = "3"                   # Human Epidermal Growth Factor 2 test
CK14_TEST = "4"                   # Cytokeratin 14 test
CK56_TEST = "5"                   # Cytokeratin 56 test

# cancer diagnoses stored in named tuple
CANCER_TYPES = ['bc1', 'bc2', 'oc', 'prc', 'pac']
CancerDiagnoses = collections.namedtuple('CancerDiagnoses', CANCER_TYPES)

# genetic tests and results stored in named tuple
GENETIC_TESTS = [gene.lower() for gene in settings.BC_MODEL['GENES']]
GeneticTests = collections.namedtuple('GeneticTests', GENETIC_TESTS)

# pathology tests stored in named tuple
PATHOLOGY_TESTS = ['er', 'pr', 'her2', 'ck14', 'ck56']
PathologyTests = collections.namedtuple('PathologyTests', PATHOLOGY_TESTS)


class PathologyTest(object):

    def __init__(self, test_type, result="0", description="pathology test"):
        self.test_type = test_type
        self.description = description
        self.result = result

    @classmethod
    def factory_default(cls):
        return PathologyTests(
            er=PathologyTest(ESTROGEN_RECEPTOR_TEST, "0", "Estrogen Receptor"),
            pr=PathologyTest(PROGESTROGEN_RECEPTOR_TEST, "0", "Progestrogen Receptor"),
            her2=PathologyTest(HER2_TEST, "0", "Human Epidermal Growth Factor 2"),
            ck14=PathologyTest(CK14_TEST, "0", "Cytokeratin 14"),
            ck56=PathologyTest(CK56_TEST, "0", "Cytokeratin 56"))

    @classmethod
    def write(cls, tests):
        """
        @type tests: PathologyTests
        @keyword tests: pathology tests
        """
        return ("%1s %1s %1s %1s %1s " %
                (tests.er.get_pathology_data(),
                 tests.pr.get_pathology_data(),
                 tests.her2.get_pathology_data(),
                 tests.ck14.get_pathology_data(),
                 tests.ck56.get_pathology_data()))

    @classmethod
    def validate(cls, person):
        """ Validate pathology data. """
        warnings = []
        tests = person.pathology
        for t in tests:
            # Check that the pathology results are correctly set (0, N, P)
            if not REGEX_PATHOLOGY_STATUS.match(t.result):
                raise PathologyError("Family member '" + person.pid + "' has been assigned an invalid " + t.test_type +
                                     " status. It must be 'N' for negative, 'P' for positive, or '0' for unknown.")
            # Check that pathology test results are only provided for family members with a first breast cancer
            if t.result != '0' and person.cancers.diagnoses.bc1.age == '0':
                raise PathologyError("Family member '" + person.pid + "' has not developed breast cancer but has " +
                                     "been assigned a breast cancer pathology test result (" + t.test_type + "). " +
                                     "Pathology test results can only be assigned to family members who have " +
                                     "developed breast cancer.")

        # if the individual has had breast cancer
        if person.cancers.diagnoses.bc1.age != "0":
            rules = "Please note the following rules for breast cancer pathology data: " \
                    "(1) if an individual's ER status is unspecified, no pathology information for that individual " \
                    "will be taken into account in the calculation; " \
                    "(2) if a breast cancer is ER positive, no other pathology information for that individual will " \
                    "be taken into account in the calculation; " \
                    "(3) if a breast cancer is ER negative, information on PR and HER2 for that individual will only " \
                    "be taken into account in the calculation if both PR and HER2 are specified; and " \
                    "(4) an individual's CK14 and CK5/6 status will only be taken into account in the calculation if " \
                    "both CK14 and CK5/6 are specified and the breast cancer is triple negative (ER negative, PR " \
                    "negative and HER2 negative). "
            # If ER is unspecified but another pathology parameter has been specified,
            # report that no pathology data will be used
            if(tests.er.result == "0" and (tests.pr.result != "0" or tests.her2.result != "0" or
                                           tests.ck14.result != "0" or tests.ck56.result != "0")):
                warnings.append(
                    "Incomplete data record in the pedigree: family member '" + person.pid + "' has an unspecified ER "
                    "status, but another pathology parameter (PR, HER2, CK14 or CK5/6) has been specified. " + rules +
                    "As a result, this individual's pathology information will not be taken into account in this case.")

            # If ER negative and PR status is specified but HER2 status is unspecified (or vice versa) report a warning
            if(tests.er.result == "N" and
               (tests.pr.result != "0" and tests.her2.result == "0") or
               (tests.pr.result == "0" and tests.her2.result != "0")):
                warnings.append(
                    "Incomplete data record in the pedigree: family member '" + person.pid + "' has a breast cancer "
                    "pathology where PR status is specified but HER2 status is unspecified (or vice versa). " + rules +
                    "As a result, PR and HER2 status will not be taken into account in this case.")

            # If either CK14 or CK5/6 has been specified (one without the other) generate a warning
            if((tests.ck14.result != "0" and tests.ck56.result == "0") or
               (tests.ck14.result == "0" and tests.ck56.result != "0")):
                warnings.append(
                    "Incomplete data record in the pedigree: family member '" + person.pid + "' has a breast cancer "
                    "pathology where only CK14 or CK5/6 status has been specified. " + rules + "As a result, CK14 and "
                    "CK5/6 status will not be taken into account in this case.")

            # If not Triple Negative but CK14 and CK5/6 are specified generate a warning
            if((tests.er.result != "N" or tests.pr.result != "N" or tests.her2.result != "N") and
               (tests.ck14.result != "0" and tests.ck56.result != "0")):
                warnings.append(
                    "Incomplete data record in your pedigree: family member '" + person.pid + "' has a breast cancer "
                    "pathology where CK14 or CK5/6 status is specified but the breast cancer pathology is not triple "
                    "negative (ER negative, PR negative and HER2 negative). " + rules + "As a result, CK14 and CK5/6 "
                    "status will not be taken into account in this case.")

            # If ER positive, and PR or HER2 or CK14 or CK5/6 specified generate a warning
            if(tests.er.result == "P" and (tests.pr.result != "0" or tests.her2.result != "0" or
                                           tests.ck14.result != "0" or tests.ck56.result != "0")):
                warnings.append(
                    "Incomplete data record in your pedigree: family member '" + person.pid + "' has a breast cancer "
                    "pathology that is ER positive, where an additional pathology parameter (PR, HER2, CK14 or CK5/6) "
                    "has been specified. " + rules + "As a result, only ER positive status will be taken into account "
                    "in this case.")

        return warnings

    def get_pathology_data(self):
        """
        Get pathology test result: '0' for untested, 'N' for negative or 'P' for positive
        @return: ER, PROG, HER2 status: neg = 1, pos = 0, unknown = 9
                 CK14, CK56 status: pos = 1, neg = 0, unknown = 9
        """
        if((not REGEX_PATHOLOGY_TEST_OPTION.match(self.test_type)) or
           (not REGEX_PATHOLOGY_STATUS.match(self.result))):
            raise PedigreeFileError(
                    "Invalid BOADICEA import format pathology test option has unexpected characters.")

        if((self.test_type == ESTROGEN_RECEPTOR_TEST) or
           (self.test_type == PROGESTROGEN_RECEPTOR_TEST) or
           (self.test_type == HER2_TEST)):
            if self.result == "0":
                return "9"    # individual is untested
            elif self.result == "N":
                return "1"    # individual has tested negative
            elif self.result == "P":
                return "0"    # individual has tested positive
            else:
                raise PedigreeFileError(
                    "Program string has unexpected value")
        elif((self.test_type == CK14_TEST) or (self.test_type == CK56_TEST)):
            if self.result == "0":
                return "9"    # individual is untested
            elif self.result == "N":
                return "0"    # individual has tested negative
            elif self.result == "P":
                return "1"    # individual has tested positive
            else:
                raise PedigreeFileError(
                    "Program string has unexpected value")

        raise PedigreeFileError(
                    "Program string has unexpected value")


class GeneticTest(object):

    def __init__(self, test_type="0", result="0"):
        """
        Genetic test.
        @keyword test_type: genetic test type ('0', 'S' or 'T')
        '0' for untested, 'S' for mutation search, or 'T' for direct gene test
        @keyword result: genetic test result ('0', 'N' or 'P')
        '0' for untested, 'N' for negative or 'P' for postive
        """
        self.test_type = test_type
        self.result = result

    @classmethod
    def validate(cls, person):
        """ Validate genetic test data. """
        gtests = person.gtests
        for t in gtests:
            # Check that the genetic test type is valid
            if not REGEX_GENETIC_TEST_TYPE.match(t.test_type):
                raise GeneticTestError("Family member '" + person.pid + "' has been assigned an invalid "
                                       "genetic test type. It must be specified with '0' for untested, "
                                       "'S' for mutation search or 'T' for direct gene test.")
            # Check that the mutation status is valid
            if not REGEX_BOADICEA_FORMAT_4_GENETIC_TEST_RESULT.match(t.result):
                raise GeneticTestError("Family member '" + person.pid + "' has been assigned an invalid "
                                       "genetic test result. Genetic test results must be '0' for untested, "
                                       "'N' for no mutation, 'P' mutation detected.")
            # If tested, check that there us a test result
            if t.test_type != "0" and t.result == '0':
                raise GeneticTestError("Family member '" + person.pid + "' has had a genetic test but the "
                                       "corresponding test result has not been specified.")
            # If there is a genetic test result check the test type is specified
            if t.test_type == "0" and t.result != '0':
                raise GeneticTestError("Family member '" + person.pid + "' has been assigned a genetic test " +
                                       "result, but the corresponding genetic test type has not been specified.")

    @classmethod
    def compareTestResults(cls, person1, person2):
        """
        Compare genetic tests of two individuals.
        @return True only if results are identical (if both have been tested) for all genes else return False
        """
        gtests1 = person1.gtests
        gtests2 = person2.gtests
        for i in range(len(GENETIC_TESTS)):
            if(REGEX_GENETIC_TEST_TYPE_IS_TESTED.match(gtests1[i].test_type) and
               REGEX_GENETIC_TEST_TYPE_IS_TESTED.match(gtests2[i].test_type)):
                if gtests1[i].result != gtests2[i].result:
                    return False
        return True

    def get_genetic_test_data(self):
        """
        Get genetic test data in 'Mendel'/Fortran input format.
        @return: 1 mutation test, mutation detected
                 2 gene test, no mutation detected
                 3 gene test, mutation detected
                 4 untested
        """
        ttype = self.test_type
        result = self.result
        if((ttype == '0') and (result == '0')):
            return 4  # untested

        # Invalid genetic test type and genetic test result
        if(((ttype == 'S') and (result == '0')) or   # (S,0) 'S' for mutation search
           ((ttype == 'T') and (result == '0')) or   # (T,0) 'T' for direct gene test
           ((ttype == '0') and (result == 'N')) or   # (0,N) untested and test result is 'N' for -ve
           ((ttype == '0') and (result == 'P'))):    # (0,P) unknown and genetic test result is 'P' for +ve
            raise GeneticTestError(
                        "Invalid BOADICEA format four genetic test summary.")

        if((ttype == 'S') and (result == 'N')):
            return '0'   # 0 mutation search test, no mutation detected

        if((ttype == 'T') and (result == 'N')):
            return '2'   # 2 gene test, no mutation detected

        if((ttype == 'S') and (result == 'P')):
            return '1'   # 1 mutation test, mutation detected

        if((ttype == 'T') and (result == 'P')):
            return '3'   # 3 gene test, mutation detected
        raise GeneticTestError(
                        "Invalid BOADICEA format four genetic test summary.")


class Cancer(object):
    """
    Basic object for cancer.
    """

    def __init__(self, age="0"):
        self.age = age


class Cancers(object):
    """
    Store diagnosis for each cancer and age of last follow up.
    """

    def __init__(self, diagnoses=CancerDiagnoses._make([Cancer("0") for _i in range(len(CANCER_TYPES))])):
        """
        @keyword diagnoses:  cancer diagnoses
        """
        self.diagnoses = diagnoses

    @classmethod
    def validate(cls, person):
        from bws import pedigree
        diagnoses = person.cancers.diagnoses
        REGEX_AGE = pedigree.REGEX_AGE
        for idx, ctype in enumerate(CANCER_TYPES):
            # Check that the age at cancer diagnosis is an unsigned integer or set to 'AU'
            # and is within range i.e. 0-110 (zero for unaffected)
            if((not REGEX_AGE.match(diagnoses[idx].age) and diagnoses[idx].age != 'AU') or
               (REGEX_AGE.match(diagnoses[idx].age) and int(diagnoses[idx].age) > settings.MAX_AGE)):
                raise CancerError("Family member '" + person.pid + "' has an age at cancer diagnosis (" + ctype +
                                  ")specified as '"+diagnoses[idx].age+"'. Age at cancer diagnosis " +
                                  "must be set to '0' for unaffected, 'AU' for affected at unknown age, or " +
                                  "specified with an integer in the range 1-"+str(settings.MAX_AGE)+".")

            # Check that the age at last follow up is greater or equal to that of all cancer diagnoses
            if(diagnoses[idx].age != 'AU' and
               diagnoses[idx].age != '0' and
               int(person.age) < int(diagnoses[idx].age)):
                raise CancerError("Family member '" + person.pid + "' has been assigned an age at cancer " +
                                  "diagnosis that exceeds age at last follow up. An age at cancer " +
                                  "diagnosis must not exceed an age at last follow up.")

        # Check that individuals who have cancer have both a valid age at last follow up and year of birth
        if person.cancers.is_cancer_diagnosed():
            # if person.cancers.alf == '0':
            #    raise CancerError("Family member '" + person.pid + "' has been diagnosed with cancer but " +
            #                      "has no age at last follow up specified. All family members with " +
            #                      "cancer must have a valid age at last follow up. If an affected " +
            #                      "family member's age is unknown, it is always better to provide " +
            #                      "some estimate of it so that risks are not underestimated.")
            if person.yob == '0':
                raise CancerError("Family member '" + person.pid + "' has been diagnosed with cancer but " +
                                  "has no year of birth specified. All family members with cancer must " +
                                  "have a valid year of birth. If an affected family member's year of " +
                                  "birth is unknown, it is always better to provide some estimate of " +
                                  "it so that risks are not underestimated.")
        # Check that males don't have an ovarian cancer diagnosis
        if(person.sex() == 'M' and diagnoses.oc.age != '0'):
            raise CancerError("Family member '" + person.pid + "' is male but has been assigned an " +
                              "ovarian cancer diagnosis.")

        # Check that females don't have a prostate cancer diagnosis
        if(person.sex() == 'F' and diagnoses.prc.age != '0'):
            raise CancerError("Family member '" + person.pid + "' is female but has been assigned an " +
                              "prostate cancer diagnosis.")

        # Check that the age of a second breast cancer exceeds that of the first.
        if(REGEX_AGE.match(diagnoses.bc2.age) and diagnoses.bc2.age != '0'):
            if diagnoses.bc1.age == '0':
                raise CancerError("Family member '" + person.pid + "' has had contralateral breast cancer, " +
                                  "but the age at diagnosis of the first breast cancer is missing.")
            elif(REGEX_AGE.match(diagnoses.bc1.age) and
                 int(diagnoses.bc1.age) > int(diagnoses.bc2.age)):
                raise CancerError("Family member '" + person.pid + "' has had contralateral breast cancer, " +
                                  "but the age at diagnosis of the first breast cancer exceeds that " +
                                  "of the second breast cancer.")

        # Check that a 2BC set to affected unknown (AU) is accompanied by a 1BC
        if(diagnoses.bc2.age == 'AU' and diagnoses.bc1.age == '0'):
            raise CancerError("Family member '" + person.pid + "' has had contralateral breast cancer, " +
                              "but the age at diagnosis of the first breast cancer is missing.")

    def write(self):
        """
        Returns a string of cancer ages used in the input pedigree file for fortran.
        """
        d = self.diagnoses
        ages = ["%3s " % c.age for c in d]
        return "".join(ages)

    def is_cancer_diagnosed(self):
        """
        Has a cancer been diagnosed.
        """
        d = self.diagnoses
        for c in d:
            if(c.age != "0"):
                return True
        return False
