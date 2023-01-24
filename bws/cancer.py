"""
Cancer, pathology and genetic testing

© 2022 Cambridge University
SPDX-FileCopyrightText: 2022 Cambridge University
SPDX-License-Identifier: GPL-3.0-or-later
"""
from collections import namedtuple
import re
import bws.consts as consts
from django.conf import settings
from django.utils.translation import gettext_lazy as _
from bws.exceptions import GeneticTestError, PathologyError, CancerError


class Genes():

    @staticmethod
    def get_unique_oc_genes():
        ''' Return genes unique to ovarian model. '''
        return list(set(settings.OC_MODEL['GENES']) - set(settings.BC_MODEL['GENES']))

    @staticmethod
    def get_all_model_genes():
        ''' Return genes for breast and ovarian model. '''
        return settings.BC_MODEL['GENES'] + settings.OC_MODEL['GENES'][4:5]


# BC pathology tests stored in named tuple
PATHOLOGY_TESTS = ['er', 'pr', 'her2', 'ck14', 'ck56']
PathologyTests = namedtuple('PathologyTests', PATHOLOGY_TESTS)


class PathologyTest(object):
    """
    Stores BC pathology test type and result.
    """
    # test types
    ESTROGEN_RECEPTOR_TEST = "1"      # Estrogen Receptor test
    PROGESTROGEN_RECEPTOR_TEST = "2"  # Progestrogen Receptor test
    HER2_TEST = "3"                   # Human Epidermal Growth Factor 2 test
    CK14_TEST = "4"                   # Cytokeratin 14 test
    CK56_TEST = "5"                   # Cytokeratin 56 test

    REGEX_PATHOLOGY_TEST_OPTION = re.compile("^([1-5])$")
    REGEX_PATHOLOGY_STATUS = re.compile("^[0NP]$")

    def __init__(self, test_type, result="0", description="pathology test"):
        self.test_type = test_type
        self.description = description
        self.result = result

    @classmethod
    def factory_default(cls):
        return PathologyTests(
            er=PathologyTest(cls.ESTROGEN_RECEPTOR_TEST, "0", "Estrogen Receptor"),
            pr=PathologyTest(cls.PROGESTROGEN_RECEPTOR_TEST, "0", "Progestrogen Receptor"),
            her2=PathologyTest(cls.HER2_TEST, "0", "HER2"),
            ck14=PathologyTest(cls.CK14_TEST, "0", "Cytokeratin 14"),
            ck56=PathologyTest(cls.CK56_TEST, "0", "Cytokeratin 5/6"))

    @classmethod
    def write(cls, tests):
        """
        Returns a string of pathology data used in the input pedigree file for fortran.
        @type tests: PathologyTests
        @keyword tests: pathology tests
        """
        pcode = -1                  # unknown or NA
        if tests.er.result == "P":
            pcode = 0               # ER+
        if tests.er.result == "N":
            if((tests.pr.result == "0" or tests.her2.result == "0") and
               (tests.pr.result != "P" and tests.her2.result != "P")):
                pcode = 1           # ER-, triple negative (TN) unknown
            elif tests.pr.result == "P" or tests.her2.result == "P":
                pcode = 2           # ER-, not TN
            elif tests.pr.result == "N" and tests.her2.result == "N":
                if tests.ck14.result == "0" or tests.ck56.result == "0":
                    pcode = 3       # TN and CK14 and/or CK56 unknown
                elif tests.ck14.result == "N" and tests.ck56.result == "N":
                    pcode = 4       # TN and CK14- and CK56-
                elif tests.ck14.result == "P" and tests.ck56.result == "P":
                    pcode = 6       # TN and CK14+ and CK56+
                elif tests.ck14.result == "P" and tests.ck56.result == "N":
                    pcode = 5   # TN and either but not both CK14+ or CK56+
                elif tests.ck14.result == "N" and tests.ck56.result == "P":
                    pcode = 5   # TN and either but not both CK14+ or CK56+
        return ("%2s " % pcode)

    @classmethod
    def validate(cls, person):
        """ Validate pathology data. """
        warnings = []
        tests = person.pathology
        for t in tests:
            # Check that the pathology results are correctly set (0, N, P)
            if not PathologyTest.REGEX_PATHOLOGY_STATUS.match(t.result):
                raise PathologyError("Family member '" + person.pid + "' has been assigned an invalid " + t.test_type +
                                     " status. It must be 'N' for negative, 'P' for positive, or '0' for unknown.",
                                     person.famid)
            # Check that pathology test results are only provided for family members with a first breast cancer
            if t.result != '0' and person.cancers.diagnoses.bc1.age == '-1':
                raise PathologyError("Family member '" + person.pid + "' has not developed breast cancer but has " +
                                     "been assigned a breast cancer pathology test result (" + t.test_type + "). " +
                                     "Pathology test results can only be assigned to family members who have " +
                                     "developed breast cancer.", person.famid)

        # if the individual has had breast cancer
        if person.cancers.diagnoses.bc1.age != "-1":
            rules = _("Please note the following rules for breast cancer pathology data: "
                      "(1) if an individual's ER status is unspecified, no pathology information for that individual "
                      "will be taken into account in the calculation; "
                      "(2) if a breast cancer is ER positive, no other pathology information for that individual will "
                      "be taken into account in the calculation; "
                      "(3) if a breast cancer is ER negative, information on PR and HER2 for that individual will only "
                      "be taken into account in the calculation if both PR and HER2 are specified; and "
                      "(4) an individual's CK14 and CK5/6 status will only be taken into account in the calculation if "
                      "both CK14 and CK5/6 are specified and the breast cancer is triple negative (ER negative, PR "
                      "negative and HER2 negative). ")
            # If ER is unspecified but another pathology parameter has been specified,
            # report that no pathology data will be used
            if(tests.er.result == "0" and (tests.pr.result != "0" or tests.her2.result != "0" or
                                           tests.ck14.result != "0" or tests.ck56.result != "0")):
                warnings.append(
                    _("Incomplete data record in the pedigree: family member \"%(id)s\" has an unspecified ER status, "
                      "but another pathology parameter (PR, HER2, CK14 or CK5/6) has been specified. %(rules)s As a "
                      "result, this individual's pathology information will not be taken into account in this case.")
                    % {'id': person.pid, 'rules': rules})

            # If ER negative and PR status is specified but HER2 status is unspecified (or vice versa) report a warning
            if(tests.er.result == "N" and
               (tests.pr.result != "0" and tests.her2.result == "0") or
               (tests.pr.result == "0" and tests.her2.result != "0")):
                warnings.append(
                    _("Incomplete data record in the pedigree: family member \"%(id)s\" has a breast cancer "
                      "pathology where PR status is specified but HER2 status is unspecified (or vice versa). %(rules)s"
                      "As a result, PR and HER2 status will not be taken into account in this case.")
                    % {'id': person.pid, 'rules': rules})

            # If either CK14 or CK5/6 has been specified (one without the other) generate a warning
            if((tests.ck14.result != "0" and tests.ck56.result == "0") or
               (tests.ck14.result == "0" and tests.ck56.result != "0")):
                warnings.append(
                    _("Incomplete data record in the pedigree: family member \"%(id)s\" has a breast cancer "
                      "pathology where only CK14 or CK5/6 status has been specified. %(rules)s As a result, CK14 and "
                      "CK5/6 status will not be taken into account in this case.")
                    % {'id': person.pid, 'rules': rules})

            # If not Triple Negative but CK14 and CK5/6 are specified generate a warning
            if((tests.er.result != "N" or tests.pr.result != "N" or tests.her2.result != "N") and
               (tests.ck14.result != "0" and tests.ck56.result != "0")):
                warnings.append(
                    _("Incomplete data record in your pedigree: family member \"%(id)s\" has a breast cancer "
                      "pathology where CK14 or CK5/6 status is specified but the breast cancer pathology is not triple "
                      "negative (ER negative, PR negative and HER2 negative). %(rules)s As a result, CK14 and CK5/6 "
                      "status will not be taken into account in this case.")
                    % {'id': person.pid, 'rules': rules})

            # If ER positive, and PR or HER2 or CK14 or CK5/6 specified generate a warning
            if(tests.er.result == "P" and (tests.pr.result != "0" or tests.her2.result != "0" or
                                           tests.ck14.result != "0" or tests.ck56.result != "0")):
                warnings.append(
                    _("Incomplete data record in your pedigree: family member \"%(id)s\" has a breast cancer pathology"
                      "that is ER positive, where an additional pathology parameter (PR, HER2, CK14 or CK5/6) "
                      "has been specified. %(rules)s As a result, only ER positive status will be taken into account "
                      "in this case.")
                    % {'id': person.pid, 'rules': rules})

        return warnings


class GeneticTest(object):
    """
    Store genetic test type and result.
    """

    REGEX_GENETIC_TEST_TYPE = re.compile("^[0ST]$")
    REGEX_BOADICEA_FORMAT_4_GENETIC_TEST_RESULT = re.compile("^[0NP]$")
    REGEX_GENETIC_TEST_TYPE_IS_TESTED = re.compile("^[ST]$")

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
            if not GeneticTest.REGEX_GENETIC_TEST_TYPE.match(t.test_type):
                raise GeneticTestError(_("Family member \"%(id)s\" has been assigned an invalid "
                                         "genetic test type. It must be specified with '0' for untested, "
                                         "'S' for mutation search or 'T' for direct gene test.")
                                       % {'id': person.pid}, person.famid)
            # Check that the mutation status is valid
            if not GeneticTest.REGEX_BOADICEA_FORMAT_4_GENETIC_TEST_RESULT.match(t.result):
                raise GeneticTestError(_("Family member \"%(id)s\" has been assigned an invalid "
                                         "genetic test result. Genetic test results must be '0' for untested, "
                                         "'N' for no mutation, 'P' mutation detected.")
                                       % {'id': person.pid}, person.famid)
            # If tested, check that there us a test result
            if t.test_type != "0" and t.result == '0':
                raise GeneticTestError(_("Family member \"%(id)s\" has had a genetic test but the "
                                         "corresponding test result has not been specified.")
                                       % {'id': person.pid}, person.famid)
            # If there is a genetic test result check the test type is specified
            if t.test_type == "0" and t.result != '0':
                raise GeneticTestError(_("Family member \"%(id)s\" has been assigned a genetic test " +
                                       "result, but the corresponding genetic test type has not been specified.")
                                       % {'id': person.pid}, person.famid)

    @classmethod
    def compareTestResults(cls, person1, person2):
        """
        Compare genetic tests of two individuals.
        @return True only if results are identical (if both have been tested) for all genes else return False
        """
        gtests1 = person1.gtests
        gtests2 = person2.gtests
        if len(gtests1) != len(gtests2):
            return False

        for idx, t1 in enumerate(gtests1):
            t2 = gtests2[idx]
            if(GeneticTest.REGEX_GENETIC_TEST_TYPE_IS_TESTED.match(t1.test_type) and
               GeneticTest.REGEX_GENETIC_TEST_TYPE_IS_TESTED.match(t2.test_type)):
                if t1.result != t2.result:
                    return False
        return True

    def get_genetic_test_data(self):
        """
        Get genetic test data in 'Mendel'/Fortran input format.
        @return: 0 mutation test, no mutation detected
                 1 mutation test, mutation detected
                 2 gene test, no mutation detected
                 3 gene test, mutation detected
                -1 untested, unknown or not applicable
        """
        ttype = self.test_type
        result = self.result
        if((ttype == '0') and (result == '0')):
            return -1  # untested

        # Invalid genetic test type and genetic test result
        if(((ttype == 'S') and (result == '0')) or   # (S,0) 'S' for mutation search
           ((ttype == 'T') and (result == '0')) or   # (T,0) 'T' for direct gene test
           ((ttype == '0') and (result == 'N')) or   # (0,N) untested and test result is 'N' for -ve
           ((ttype == '0') and (result == 'P'))):    # (0,P) unknown and genetic test result is 'P' for +ve
            raise GeneticTestError("Invalid BOADICEA format for genetic test summary.")

        if((ttype == 'S') and (result == 'N')):
            return '0'   # 0 mutation search test, no mutation detected

        if((ttype == 'T') and (result == 'N')):
            return '2'   # 2 gene test, no mutation detected

        if((ttype == 'S') and (result == 'P')):
            return '1'   # 1 mutation test, mutation detected

        if((ttype == 'T') and (result == 'P')):
            return '3'   # 3 gene test, mutation detected
        raise GeneticTestError("Invalid BOADICEA format for genetic test summary.")


class GeneticTestsMixin():
    __slots__ = ()   # set to empty tuple prevents creation of instance dictionary and keeps memory low

    @classmethod
    def default_factory(cls):
        """ Set all genetic tests to untested """
        return cls._make([GeneticTest() for _i in range(len(cls._fields))])

    @classmethod
    def factory(cls, geneticTestList):
        return cls._make(geneticTestList)


class BWSGeneticTests(namedtuple('BWSGeneticTests', ' '.join([gene.lower() for gene in settings.BC_MODEL['GENES']])),
                      GeneticTestsMixin):
    """
    BWS Genetic tests and results stored in named tuple of gene names
    """


class CanRiskGeneticTests(namedtuple('CanRiskGeneticTests',
                                     ' '.join([gene.lower() for gene in Genes.get_all_model_genes()])
                                     ),
                          GeneticTestsMixin):
    """
    CanRisk Genetic tests and results stored in named tuple of gene names
    """


class Cancer(object):
    """
    Basic object for cancer.
    """
    def __init__(self, age="-1"):
        self.age = age


class Cancers():
    """
    Store diagnosis for each cancer and age of last follow up.
    """

    def __init__(self, **kwargs):
        """
        @keyword kwargs: list of Cancer objects
        """
        cancer_types = Cancers.get_cancers()
        for ctype in cancer_types:
            if ctype not in kwargs:
                kwargs[ctype] = Cancer()

        # cancer diagnoses stored in named tuple
        CancerDiagnoses = namedtuple('CancerDiagnoses', cancer_types)
        self.diagnoses = CancerDiagnoses(**kwargs)

    @classmethod
    def validate(cls, person):
        """
        Validate a person's cancer types and diagnoses ages.
        """
        cancer_types = Cancers.get_cancers()
        diagnoses = person.cancers.diagnoses
        REGEX_AGE = consts.REGEX_AGE
        for idx, ctype in enumerate(cancer_types):
            diagnoses_age = diagnoses[idx].age
            # Check that the age at cancer diagnosis is an unsigned integer or set to 'AU'
            # and is within range i.e. 0-110 (-1 for unaffected)
            if((not REGEX_AGE.match(diagnoses_age) and diagnoses_age != 'AU' and diagnoses_age != '-1') or
               (REGEX_AGE.match(diagnoses_age) and int(diagnoses_age) > settings.MAX_AGE)):
                    raise CancerError(_("Family member \"%(id)s\" has an age at cancer diagnosis (%(ctype)s) "
                                        "specified as %(dage)s. Age at cancer diagnosis " +
                                        "must be set to '0' for unaffected, 'AU' for affected at unknown age, or " +
                                        "specified with an integer in the range 1-%(max_age)s.")
                                      % {'id': person.pid, 'ctype': ctype, 'dage': diagnoses_age,
                                         'max_age': settings.MAX_AGE}, person.famid)

            # Check that the age at last follow up is greater or equal to that of all cancer diagnoses
            if(diagnoses_age != 'AU' and diagnoses_age != '-1' and int(person.age) < int(diagnoses_age)):
                raise CancerError(_("Family member \"%(id)s\" has been assigned an age at cancer " +
                                    "diagnosis that exceeds age at last follow up. An age at cancer " +
                                    "diagnosis must not exceed an age at last follow up.")
                                  % {'id': person.pid}, person.famid)

            # Check that males don't have an ovarian cancer diagnosis
            if ctype == 'oc' and person.sex() == 'M' and diagnoses_age != '-1':
                raise CancerError(_("Family member \"%(id)s\" is male but has been assigned an " +
                                    "ovarian cancer diagnosis.")
                                  % {'id': person.pid}, person.famid)

            # Check that females don't have a prostate cancer diagnosis
            if ctype == 'prc' and person.sex() == 'F' and diagnoses_age != '-1':
                raise CancerError(_("Family member \"%(id)s\" is female but has been assigned an " +
                                    "prostate cancer diagnosis.")
                                  % {'id': person.pid}, person.famid)

        # Check that individuals who have cancer have both a valid age at last follow up and year of birth
        cancers = person.cancers
        if cancers.is_cancer_diagnosed():
            if person.yob == '0':
                raise CancerError(_("Family member \"%(id)s\" has been diagnosed with cancer but " +
                                    "has no year of birth specified. All family members with cancer must " +
                                    "have a valid year of birth. If an affected family member's year of " +
                                    "birth is unknown, it is always better to provide some estimate of " +
                                    "it so that risks are not underestimated.")
                                  % {'id': person.pid}, person.famid)

        # Check that the age of a second breast cancer exceeds that of the first.
        bc1 = getattr(cancers.diagnoses, "bc1", None)
        bc2 = getattr(cancers.diagnoses, "bc2", None)
        if bc1 is not None and bc2 is not None:
            if(REGEX_AGE.match(bc2.age) and bc2.age != '-1'):
                if bc1.age == '-1':
                    raise CancerError(_("Family member \"%(id)s\" has had contralateral breast cancer, " +
                                        "but the age at diagnosis of the first breast cancer is missing.")
                                      % {'id': person.pid}, person.famid)
                elif(REGEX_AGE.match(bc1.age) and int(bc1.age) > int(bc2.age)):
                    raise CancerError(_("Family member \"%(id)s\" has had contralateral breast cancer, " +
                                        "but the age at diagnosis of the first breast cancer exceeds that " +
                                        "of the second breast cancer.")
                                      % {'id': person.pid}, person.famid)

            # Check that a 2BC set to affected unknown (AU) is accompanied by a 1BC
            if(bc2.age == 'AU' and bc1.age == '-1'):
                raise CancerError(_("Family member \"%(id)s\" has had contralateral breast cancer, " +
                                    "but the age at diagnosis of the first breast cancer is missing.")
                                  % {'id': person.pid}, person.famid)

    def write(self, cancers=None, age=-1):
        """
        Returns a string of cancer ages used in the input pedigree file for fortran.
        Affected Unknown (AU) are set to the age of last followup to avoid underestimation of risk.
        """
        d = self.diagnoses
        ages = ["%3s " % (getattr(d, c).age if getattr(d, c).age != 'AU' else age) for c in cancers]
        return "".join(ages)

    def is_cancer_diagnosed(self):
        """
        Has a cancer been diagnosed.
        """
        d = self.diagnoses
        for c in d:
            if c.age != "-1":
                return True
        return False

    @classmethod
    def get_cancers(cls):
        """ Get a list of the cancer types stored for CanRisk. """
        return ['bc1', 'bc2', 'oc', 'prc', 'pac']
