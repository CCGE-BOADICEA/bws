"""
Person in the Pedigree

© 2023 University of Cambridge
SPDX-FileCopyrightText: 2023 University of Cambridge
SPDX-License-Identifier: GPL-3.0-or-later
"""

from datetime import date
from django.conf import settings

from bws.cancer import Cancer, GeneticTest, PathologyTests, PathologyTest, Cancers, \
    BWSGeneticTests, CanRiskGeneticTests, Genes
import bws.consts as consts
from bws.exceptions import PedigreeError, PersonError
import bws.pedigree as pedigree


class Person(object):
    """ Person class. """

    def __init__(self, famid, name, pid, fathid, mothid, target="0", dead="0", age="0", yob="0", ashkn="0", mztwin="0",
                 cancers=Cancers(),
                 gtests=BWSGeneticTests.default_factory(),
                 pathology=PathologyTest.factory_default()):
        """
        @type  famid: str
        @param famid: family/pedigree ID
        @type  pid: str
        @param pid: person ID
        @type  fathid: str
        @keyword fathid: father ID
        @type  mothid: str
        @keyword mothid: mother ID
        @type  target: str
        @keyword target: subject of risk calculation
        @type  dead: str
        @keyword dead: alive specified as '0' and dead as '1'
        @type  age: str
        @keyword age: age at last follow up or age at death
        @type  yob: str
        @keyword yob: year of birth
        @type ashkn: str
        @keyword ashkn: Ashkenazi origin parameter: '1' for Ashkenazi origin else '0'
        @type  mztwin: str
        @keyword mztwin: monozygotic (identical) twin
        @type cancers: Cancers
        @keyword cancers: cancer status
        @type gtest: GeneticTests
        @keyword gtest: genetic tests
        @type pathology: PathologyResult
        @keyword pathology: pathology test results
        """
        self.famid = famid.replace("-", "")[:8]   # remove hyphen and restrict to 8 chars
        self.name = name[:8]
        self.pid = pid
        self.fathid = fathid
        self.mothid = mothid
        self.target = target
        self.dead = dead
        self.age = age
        self.yob = yob
        self.ashkn = ashkn
        self.mztwin = mztwin
        self.cancers = cancers  # cancers
        self.gtests = gtests    # genetic tests
        self.pathology = pathology

    def validate(self, pedigree):
        """ Validation check for people input.
        @param pedigree: Pedigree the person belongs to.
        """
        if(self.name == '' or
           not consts.REGEX_ALPHANUM_HYPHENS.match(self.name)):
            raise PersonError("A name '"+self.name+"' is unspecified or is not an alphanumeric string.", self.famid)

        if(len(self.pid) < settings.MIN_FAMILY_ID_STR_LENGTH or
           len(self.pid) > settings.MAX_FAMILY_ID_STR_LENGTH or
           consts.REGEX_ONLY_ZEROS.match(self.pid) or
           not consts.REGEX_ALPHANUM_HYPHENS.match(self.pid)):
            raise PersonError("An individual identifier (IndivID column) was specified as '" + self.pid +
                              ". Individual identifiers must be alphanumeric strings with a maximum of " +
                              str(settings.MAX_FAMILY_ID_STR_LENGTH)+"characters.", self.famid)

        if(len(self.fathid) < settings.MIN_FAMILY_ID_STR_LENGTH or
           len(self.fathid) > settings.MAX_FAMILY_ID_STR_LENGTH or
           not consts.REGEX_ALPHANUM_HYPHENS.match(self.fathid)):
            raise PersonError("Father identifier ('" + self.fathid + "', FathID column) has unexpected characters. "
                              "It must be alphanumeric strings with a maximum of " +
                              str(settings.MAX_FAMILY_ID_STR_LENGTH) + " characters", self.famid)

        if(len(self.mothid) < settings.MIN_FAMILY_ID_STR_LENGTH or
           len(self.mothid) > settings.MAX_FAMILY_ID_STR_LENGTH or
           not consts.REGEX_ALPHANUM_HYPHENS.match(self.mothid)):
            raise PersonError("Mother identifier ('" + self.mothid + "', MothID column) has unexpected characters. "
                              "It must be alphanumeric strings with a maximum of " +
                              str(settings.MAX_FAMILY_ID_STR_LENGTH) + " characters", self.famid)

        if(self.fathid == '0' and self.mothid != '0') or (self.fathid != '0' and self.mothid == '0'):
            raise PersonError("Family member '"+self.name+"' has only one parent specified. All family members must "
                              "have no parents specified (i.e. they must be founders) or both parents specified.",
                              self.famid)

        # check for missing parents
        if self.mothid != '0' and pedigree.get_person(self.mothid) is None:
            raise PersonError("The mother '"+self.mothid+"' of family member '" + self.pid +
                              "' is missing from the pedigree.", self.famid)
        elif self.fathid != '0' and pedigree.get_person(self.fathid) is None:
            raise PersonError("The father '"+self.fathid+"' of family member '" + self.pid +
                              "' is missing from the pedigree.", self.famid)

        # check all fathers are male
        if self.fathid != '0' and pedigree.get_person(self.fathid).sex() != 'M':
            raise PersonError("The father of family member '" + self.pid + "' is not specified as male. " +
                              "All fathers in the pedigree must have sex specified as 'M'.", self.famid)
        # check all mothers are female
        if self.mothid != '0' and pedigree.get_person(self.mothid).sex() != 'F':
            raise PersonError("The mother of family member '" + self.pid + "' is not specified as female. " +
                              "All mothers in the pedigree must have sex specified as 'F'.", self.famid)
        # check if the dead attribute has been correctly set
        if self.dead != '0' and self.dead != '1':
            raise PersonError("The family member '" + self.pid + "' has an invalid vital status " +
                              "(alive must be specified as '0', and dead specified as '1')", self.famid)
        # check that age of last follow up set to either 0 (unknown) or in range 1-110
        if not consts.REGEX_AGE.match(self.age) or int(self.age) > settings.MAX_AGE:
            raise PersonError("The age specified for family member '" + self.pid + "' has unexpected " +
                              "characters. Ages must be specified with as '0' for unknown, or in the " +
                              "range 1-" + str(settings.MAX_AGE), self.famid)

        # validate year of birth
        current_year = date.today().year
        if self.yob != "0":
            if(not consts.REGEX_YEAR_OF_BIRTH.match(self.yob) or
               int(self.yob) < settings.MIN_YEAR_OF_BIRTH or
               int(self.yob) > current_year):
                raise PersonError("The year of birth '" + self.yob + "' specified for family member '" + self.pid +
                                  "' is out of range. Years of birth must be in the range " +
                                  str(settings.MIN_YEAR_OF_BIRTH) + "-" + str(current_year))

        # Check that we have valid data values for the Ashkenazi flag
        if(not consts.REGEX_ASHKENAZI_STATUS.match(self.ashkn)):
            raise PersonError("Family member '" + self.pid + "' has been assigned an invalid Ashkenazi "
                              "origin parameter. The Ashkenazi origin parameter must be set to '1' "
                              "for Ashkenazi origin, or '0' for not Ashkenazi origin.")

        # Check max no. of siblings not exceeded
        (siblings, siblings_same_yob) = pedigree.get_siblings(self)
        if len(siblings) > settings.MAX_NUMBER_OF_SIBS_PER_NUCLEAR_FAMILY:
            raise PersonError("Family member '" + self.pid + "' exceeded the maximum number of siblings (" +
                              str(settings.MAX_NUMBER_OF_SIBS_PER_NUCLEAR_FAMILY) + ".")
        # Check has siblings with the same year of birth
        if len(siblings_same_yob) > settings.MAX_NUMBER_OF_SIBS_PER_NUCLEAR_FAMILY_WITH_SAME_YOB:
            raise PersonError("Family member '" + self.pid + "' exceeded the maximum number of siblings " +
                              "with the same year of birth exceeded (" +
                              str(settings.MAX_NUMBER_OF_SIBS_PER_NUCLEAR_FAMILY_WITH_SAME_YOB) + ")")

    @staticmethod
    def factory(ped_file_line, file_type=None):
        ''' Factory method for creating types of people given a record from
        a BOADICEA import pedigree file .
        @type  ped_file_line: str
        @param ped_file_line: Pedigree file line.
        '''
        cols = ped_file_line.split()

        famid = cols[0]
        name = cols[1]
        pid = cols[3]
        cancers = Cancers(bc1=Cancer(cols[11] if cols[11] != "0" else "-1"),
                          bc2=Cancer(cols[12] if cols[12] != "0" else "-1"),
                          oc=Cancer(cols[13] if cols[13] != "0" else "-1"),
                          prc=Cancer(cols[14] if cols[14] != "0" else "-1"),
                          pac=Cancer(cols[15] if cols[15] != "0" else "-1"))

        # use column headers to get gene test type and result
        if file_type == 'bwa':
            gtests = BWSGeneticTests.factory([GeneticTest(cols[pedigree.BwaPedigree.get_column_idx(gene+'t')],
                                                          cols[pedigree.BwaPedigree.get_column_idx(gene+'r')])
                                             if pedigree.BwaPedigree.get_column_idx(gene+'t') != -1 else GeneticTest()
                                             for gene in settings.BC_MODEL['GENES']])
            pathology = PathologyTests(
                er=PathologyTest(PathologyTest.ESTROGEN_RECEPTOR_TEST, cols[27]),
                pr=PathologyTest(PathologyTest.PROGESTROGEN_RECEPTOR_TEST, cols[28]),
                her2=PathologyTest(PathologyTest.HER2_TEST, cols[29]),
                ck14=PathologyTest(PathologyTest.CK14_TEST, cols[30]),
                ck56=PathologyTest(PathologyTest.CK56_TEST, cols[31]))
        else:
            genes = Genes.get_all_model_genes()

            def get_genetic_test(cols, gene):
                idx = pedigree.CanRiskPedigree.get_column_idx(gene, file_type)
                if idx < 0:
                    if gene == "BARD1" and file_type == "canrisk1":
                        return GeneticTest()
                    raise PedigreeError("Genetic test column for '" + gene + "not found.")
                gt = cols[idx].split(':')
                return GeneticTest(gt[0], gt[1])
            gtests = CanRiskGeneticTests.factory([get_genetic_test(cols, gene) for gene in genes])

            path = cols[pedigree.CanRiskPedigree.get_column_idx("ER:PR:HER2:CK14:CK56", file_type)].split(':')
            pathology = PathologyTests(
                er=PathologyTest(PathologyTest.ESTROGEN_RECEPTOR_TEST, path[0]),
                pr=PathologyTest(PathologyTest.PROGESTROGEN_RECEPTOR_TEST, path[1]),
                her2=PathologyTest(PathologyTest.HER2_TEST, path[2]),
                ck14=PathologyTest(PathologyTest.CK14_TEST, path[3]),
                ck56=PathologyTest(PathologyTest.CK56_TEST, path[4]))

        if cols[6] == 'M':
            return Male(
                famid, name, pid, fathid=cols[4], mothid=cols[5], target=cols[2], dead=cols[8], age=cols[9],
                yob=cols[10], ashkn=cols[16], cancers=cancers, mztwin=cols[7], gtests=gtests, pathology=pathology)
        elif cols[6] == 'F':
            return Female(
                famid, name, pid, fathid=cols[4], mothid=cols[5], target=cols[2], dead=cols[8], age=cols[9],
                yob=cols[10], ashkn=cols[16], cancers=cancers, mztwin=cols[7], gtests=gtests, pathology=pathology)
        else:
            raise PedigreeError("The sex of family member '"+name+"' is invalid. An " +
                                "individuals sex must be specified as 'M' or 'F' only.")

    def is_complete(self):
        """
        Checks if a family member without cancer is missing a year of birth or age at last
        follow up. An individual's year of birth and age at last follow up must be specified
        in order for that person to be included in a calculation. As a result, family members
        lacking this information will be excluded from the calculation.
        """
        if self.yob == '0' or self.age == '0':
            if not self.cancers.is_cancer_diagnosed():
                return False
        return True

    def is_target(self):
        return False if self.target == '0' else True


class Male(Person):
    ''' Male person. '''

    def sex(self):
        return 'M'


class Female(Person):
    ''' Female person. '''

    def sex(self):
        return 'F'
