"""
Pedigree data

© 2022 Cambridge University
SPDX-FileCopyrightText: 2022 Cambridge University
SPDX-License-Identifier: GPL-3.0-or-later
"""
import abc
import logging
import os
from random import randint

from django.conf import settings

from bws.cancer import GeneticTest, PathologyTest, BWSGeneticTests, Genes
import bws.consts as consts
from bws.exceptions import PedigreeError
from bws.person import Person, Male, Female
from bws.risk_factors.mdensity import Volpara, Stratus


logger = logging.getLogger(__name__)


class Pedigree(metaclass=abc.ABCMeta):
    """
    A pedigree object.
    """

    def __init__(self, pedigree_records=None, people=None, file_type=None,
                 bc_risk_factor_code=None, oc_risk_factor_code=None,
                 bc_prs=None, oc_prs=None, hgt=-1, mdensity=None):
        """
        @keyword pedigree_records: the pedigree records section of the BOADICEA import pedigree file.
        @keyword people: members of the pedigree.
        @keyword file_type: file type is 'bwa' or 'canrisk'.
        @keyword bc_risk_factor_code: breast cancer risk factor code
        @keyword oc_risk_factor_code: ovarian cancer risk factor code
        @keyword bc_prs: breast cancer PRS
        @keyword oc_prs: ovarian cancer PRS
        """
        self.people = []
        if pedigree_records is not None:
            self.famid = pedigree_records[0].split()[0]
            ids = []
            for record in pedigree_records:
                p = Person.factory(record, file_type=file_type)
                if p.target != '0' and p.target != '1':
                    raise PedigreeError("A value in the Target data column has been set to '" + p.target +
                                        "'. Target column parameters must be set to '0' or '1'.", p.famid)
                if p.is_target():
                    self.target = p

                if p.pid in ids:
                    raise PedigreeError("Individual ID '" + p.pid +
                                        "' appears more than once in the pedigree file.", p.famid)
                else:
                    ids.append(p.pid)
                self.people.append(p)
        if people is not None:
            self.people.extend(people)
            self.famid = self.people[0].famid

        ntarget = 0
        for person in self.people:
            if person.is_target():
                ntarget += 1

        pedigree_size = len(self.people)
        if ntarget != 1:
            raise PedigreeError("Pedigree (" + self.famid + ") has either no index or more than 1 " +
                                "index individuals. Only one target can be specified.", self.famid)
        if pedigree_size > settings.MAX_PEDIGREE_SIZE or pedigree_size < settings.MIN_BASELINE_PEDIGREE_SIZE:
            raise PedigreeError("Pedigree (" + self.famid + ") has unexpected number of family members " +
                                str(pedigree_size), self.famid)
        if file_type is not None and file_type.startswith('canrisk'):
            self.hgt = hgt
            self.mdensity = mdensity
            if bc_risk_factor_code is not None:
                self.bc_risk_factor_code = bc_risk_factor_code
            if oc_risk_factor_code is not None:
                self.oc_risk_factor_code = oc_risk_factor_code
            if bc_prs is not None:
                self.bc_prs = bc_prs
            if oc_prs is not None:
                self.oc_prs = oc_prs

    def validate(self):
        """ Validation check for pedigree input.
        @param p: Person to validate pedigree data.
        """
        if(len(self.famid) > settings.MAX_LENGTH_PEDIGREE_NUMBER_STR or
           not consts.REGEX_ALPHANUM_HYPHENS.match(self.famid) or       # must be alphanumeric plus hyphen
           consts.REGEX_ONLY_HYPHENS.match(self.famid) or               # but not just hyphens
           consts.REGEX_ONLY_ZEROS.match(self.famid)):                  # and not just zeros
            raise PedigreeError(
                "Family ID (1st data column) has been set to '" + self.famid +
                "'. Family IDs must be specified with between 1 and "+str(settings.MAX_LENGTH_PEDIGREE_NUMBER_STR) +
                " non-zero number or alphanumeric characters.", self.famid)

        unconnected = self.unconnected()
        if len(unconnected) > 0:
            raise PedigreeError("Pedigree (" + self.famid + ") family members are not physically " +
                                "connected to the target: " + str(unconnected), self.famid)

        # Check that the index's parameters are valid
        target = self.get_target()
        if target.yob == '0':
            raise PedigreeError("The target's year of birth has been set to '" + target.yob +
                                "'. This person must be assigned a valid year of birth.", target.famid)
        if target.age == '0':
            raise PedigreeError("The target's age has been set to '" + target.age +
                                "'. This person must be assigned an age.", target.famid)

        # Check that carrier probabilities / cancer risks can be computed
        carrier_probs = self.is_carrier_probs_viable(target=target)
        cancer_risks = self.is_risks_calc_viable(target=target)
        if(not carrier_probs and not cancer_risks):
            raise PedigreeError(
                "BOADICEA cannot compute mutation carrier probabilities because the target '" + target.pid +
                "' has a positive genetic test. Also BOADICEA cannot compute breast and ovarian cancer "
                "risks because the target is: (1) over " + str(settings.MAX_AGE_FOR_RISK_CALCS) +
                " years old or (2) male, or (3) an affected female who has developed contralateral "
                "breast cancer, ovarian cancer or pancreatic cancer.", target.famid)

        #
        # Check monozygotic (MZ) twin data
        twin_store = self.get_twins()

        # Check that MZ siblings are only specified as twins, no identical triplets etc
        for t in twin_store:
            twins = twin_store[t]
            if len(twins) != 2:
                raise PedigreeError(
                    "MZ twin identifier '" + str(twins[0].pid) + "' does not appear twice in the pedigree file. "
                    "Only MZ twins are permitted in the pedigree, MZ triplets or quads are not allowed.",
                    twins[0].famid)

            # Check MZ twin characters are valid
            if len(t) != 1 or t not in settings.UNIQUE_TWIN_IDS:
                raise PedigreeError("Invalid MZ twin character '" + t + "'. MZ twins must be identified using one " +
                                    "of the following ASCII characters: " + str(settings.UNIQUE_TWIN_IDS) + ".",
                                    twins[0].famid)

            # Check that monozygotic (MZ) twin data are consistent
            if(twins[0].mothid != twins[1].mothid or
               twins[0].fathid != twins[1].fathid):
                raise PedigreeError("Monozygotic (MZ) twins identified with the character '" + t + "' have different "
                                    "parents. MZ twins must have the same parents.", twins[0].famid)
            if(twins[0].yob != twins[1].yob):
                raise PedigreeError("Monozygotic (MZ) twins identified with the character '" + t + "' have different "
                                    "years of birth. MZ twins must have the same year of birth.", twins[0].famid)

            # Check that living MZ twins have the same age at last follow up
            if(twins[0].dead == '0' and twins[1].dead == '0' and twins[0].age != twins[1].age):
                raise PedigreeError("Monozygotic (MZ) twins identified with the character '" + t + "' have different "
                                    "ages. If both MZ twins are alive, they must have the same age at last follow up.",
                                    twins[0].famid)

            if twins[0].sex() != twins[1].sex():
                raise PedigreeError("Monozygotic (MZ) twins identified with the character '" + t + "' have a different "
                                    "sex. MZ twins must have the same sex.", twins[0].famid)

            # Check that the MZ twins have the same genetic status
            if not GeneticTest.compareTestResults(twins[0], twins[1]):
                raise PedigreeError("Monozygotic (MZ) twins have both had a genetic test, but the genetic test results "
                                    "for these individuals are different. Under these circumstances, the genetic test "
                                    "results must be the same.", twins[0].famid)

        # Check to ensure that the maximum number of MZ twin pairs per pedigree has not been exceeded
        if len(twin_store.keys()) > settings.MAX_NUMBER_MZ_TWIN_PAIRS:
            raise PedigreeError("Maximum number of MZ twin pairs has been exceeded. Input pedigrees must have a "
                                "maximum of " + str(settings.MAX_NUMBER_MZ_TWIN_PAIRS) + " MZ twin pairs.", self.famid)

    def add_parents(self, person, gtests=BWSGeneticTests.default_factory()):
        """
        Add parents for a given person to the pedigree.
        @param person: Person in the pedigree to add parents to
        @keyword gtests: genetic test results
        @return: father and mother
        """
        if person.fathid == "0" and person.mothid == "0":
            n = randint(1000, 9999)
            person.fathid = person.pid + str(n)
            person.mothid = person.pid + str(n+1)
        father = Male(person.famid, person.name + "father", person.fathid, "0", "0", gtests=gtests)
        mother = Female(person.famid, person.name + "mother", person.mothid, "0", "0", gtests=gtests)
        self.people.append(father)
        self.people.append(mother)
        return (father, mother)

    def get_target(self):
        """
        Get target in the pedigree.
        @return: target
        """
        for p in self.people:
            if p.is_target():
                return p
        return None

    def get_person(self, individ):
        """
        Get a person in the pedigree by their IndivID.
        @return: the requested person
        """
        for p in self.people:
            if p.pid == individ:
                return p
        return None

    def is_ashkn(self):
        """
        Does the pedigree have Ashkenazi Jewish ancestry
        @return: True if Ashkenazi Jewish ancestry
        """
        for p in self.people:
            if p.ashkn == "1":
                return True
        return False

    def get_siblings(self, person):
        """
        Get a list of the siblings of the given person.
        """
        individ = person.pid
        fathid = person.fathid
        mothid = person.mothid
        siblings = []
        siblings_same_yob = []
        if fathid == "0" or mothid == "0":
            return (siblings, siblings_same_yob)
        for p in self.people:
            if p.pid != individ:
                if p.mothid == mothid and p.fathid == fathid:
                    siblings.append(p)
                    if p.yob == person.yob:
                        siblings_same_yob.append(p)
        return (siblings, siblings_same_yob)

    def get_person_by_name(self, name):
        """
        Get a person in the pedigree by their name.
        @return: the requested person
        """
        for p in self.people:
            if p.name == name:
                return p
        return None

    def get_twins(self):
        """
        Get a dictionary of the monozygotic (MZ) twins in the pedigree
        @return: dictionary of mztwins
        """
        twin_store = {}
        for p in self.people:
            if p.mztwin != "0":
                twin = p.mztwin
                if twin in twin_store:
                    twin_store[twin].append(p)
                else:
                    twin_store[twin] = [p]
        return twin_store

    def unconnected(self):
        """
        Based on Andrew Lee's mod_pedigree.is_connected() routine.
        Determines those people connected to the proband, and identifies those individuals
        that are not connected.
        @return: return a list of individuals that aren't connected to the target
        """
        target = self.get_target()
        connected = [target.pid]
        change = True
        while change:
            change = False
            for p in self.people:
                # If the person is in the connected group, but their parents aren't, then add them
                if p.pid in connected:
                    if((p.mothid != '0') and (p.mothid not in connected)):
                        connected.append(p.mothid)
                        change = True
                    if((p.fathid != '0') and (p.fathid not in connected)):
                        connected.append(p.fathid)
                        change = True
                # If the person parents is in the connected group, but they aren't, then add them
                elif ((p.mothid in connected) or (p.fathid in connected)):
                    connected.append(p.pid)
                    change = True
        # list of the individuals not connected to the target.
        return [p.pid for p in self.people if p.pid not in connected]

    def is_risks_calc_viable(self, target=None):
        """
        Returns true if the target is female and unaffected or only one breast cancer. If
        the target is male, female with additional cancers, age at last follow up is not
        in specified range or year of birth is not valid return false.
        @return: true if risks calculation is viable
        """
        if target is None:
            target = self.get_target()
        d = target.cancers.diagnoses
        if ((target is None or
             isinstance(target, Male) or
             (d.bc2.age != "-1") or (d.pac.age != "-1") or (d.oc.age != "-1") or
             (int(target.age) > settings.MAX_AGE_FOR_RISK_CALCS) or
             (int(target.yob) < settings.MIN_YEAR_OF_BIRTH))):
            return False
        return True

    def is_carrier_probs_viable(self, target=None, genes=Genes.get_all_model_genes()):
        """
        Return true if the target does not have a positive genetic test carrier probs
        cannot be calculated.
        @return: true if carrier probability calculation is viable
        """
        if target is None:
            target = self.get_target()
        gtests = target.gtests
        for g in genes:
            t = getattr(gtests, g.lower(), GeneticTest())
            if t.result == 'P':
                return False
        return True

    def write_pedigree_file(self, risk_factor_code='0', hgt=-1, mdensity=None, prs=None, filepath="/tmp/test.ped",
                            model_settings=settings.BC_MODEL):
        """
        Write input pedigree file for fortran.
        """
        
        if mdensity is not None and (isinstance(mdensity, Volpara) or isinstance(mdensity, Stratus)):
            raise Exception("Unsupported mammographic density type")
        
        f = open(filepath, "w")
        
        mname = model_settings['NAME']
        num = "5" if mname == "BC" else "4"  # extra column for mammographic density
        print("(I3,X,A8)", file=f)

        print("(3(A7,X),2(A1,X),2(A3,X)," + str(len(model_settings['CANCERS'])+1) + "(A3,X)," +
              str(len(model_settings['GENES'])) + "(A2,X),A4,X,A2,X,A1," + num + "(X,A8))", file=f)

        print("%-3d %-8s" % (len(self.people), self.people[0].famid), file=f)

        for p in self.people:
            # IndivID FathID MothID Sex MZ Genotype, Polygene 1BC 2BC OC
            genotype = ''
            print("%-7s %-7s %-7s %-1s %-1s %3s %-3s " %
                  (p.pid,
                   p.fathid if p.fathid != "0" else '',
                   p.mothid if p.mothid != "0" else '', p.sex(),
                   p.mztwin if p.mztwin != "0" else '',
                   genotype, '   '), file=f, end="")

            print(p.cancers.write(model_settings['CANCERS'], p.age), file=f, end="")
            print("%3s " % p.age, file=f, end="")

            # Gene Tests
            gtests = p.gtests
            for g in model_settings['GENES']:
                try:
                    print("%2s " % getattr(gtests, g.lower()).get_genetic_test_data(), file=f, end="")
                except AttributeError:
                    # check if gene not in BC model
                    if mname == "OC" and isinstance(gtests, BWSGeneticTests):
                        if g in Genes.get_unique_oc_genes():
                            print("%2s " % GeneticTest().get_genetic_test_data(), file=f, end="")
                    else:
                        raise

            print("%4s " % (p.yob if p.yob != "0" else settings.MENDEL_NULL_YEAR_OF_BIRTH), file=f, end="")

            print(PathologyTest.write(p.pathology), file=f, end="")

            # ProbandStatus RiskFactor
            print("%1s %8s " % (p.target, (risk_factor_code if p.target != "0" else "00000000")),
                  file=f, end="")

            # Height
            print(("%8.4f " % hgt) if p.target != "0" else ("%8s " % "-1"), file=f, end="")

            # Mammographic density
            if mname == "BC":
                print(("%8s " % mdensity.get_pedigree_str()) if p.target != "0" and mdensity is not None else ("%8s " % "00000000"), file=f, end="")

            # PolygStanDev PolygLoad
            print("%8.5f %8.5f" % (prs.alpha if p.target != "0" and prs is not None and prs.alpha else 0,
                                   prs.zscore if p.target != "0" and prs is not None and prs.zscore else 0,),
                  file=f)

        f.close()
        return filepath

    def write_param_file(self, filepath="/tmp/params",
                         model_settings=settings.BC_MODEL,
                         mutation_freq=settings.BC_MODEL['MUTATION_FREQUENCIES']['UK'],
                         sensitivity=settings.BC_MODEL['GENETIC_TEST_SENSITIVITY'],
                         isashk=False):
        """
        Write model parameters file.
        @param filepath: path to write the model parameters file to
        @param model_settings: model settings
        @param mutation_freq: mutation frequencies
        @param sensitivity: genetic test sensitivity
        @param isashk: true if AJ
        """
        # Note: population allele frequencies are used to compute the incidence rates
        # for each genotype, from the overall population incidences
        allele_freq = "PEDIGREE_ALLELE_FRQ" if isashk else "POPULATION_ALLELE_FRQ"
        f = open(filepath, "w")
        print("&settings", file=f, end='\n\n')
        for idx, gene in enumerate(model_settings['GENES'], start=1):
            print(f"{allele_freq}( {idx} ) = {mutation_freq[gene]}", file=f)
        for idx, gene in enumerate(model_settings['GENES'], start=1):
            print(f"SCREENING_SENSITIVITIES( {idx} ) = {sensitivity[gene]}", file=f)
        print("/", file=f)
        f.close()
        return filepath

    def write_batch_file(self, pedigree_file_name, filepath="/tmp/test.bat",
                         model_settings=settings.BC_MODEL, calc_ages=None):
        """
        Write fortran input batch file.
        @param pedigree_file_name: path to fortran pedigree file
        @param filepath: path to write the batch file to
        @param model_settings: model settings
        @param calc_ages: list of ages to calculate a cancer risk at
        """
        f = open(filepath, "w")

        print("2", file=f)
        print(os.path.join(model_settings['HOME'], "Data/locus.loc"), file=f)

        target = self.get_target()
        tage = int(target.age)      # target age at last follow up

        if calc_ages is None:
            # Compute breast/ovarian cancer risks for the following years:
            #    (1) Next 5 years at one year intervals, age at last follow up +1, +2, +3, +4, +5
            #    (2) Age at last follow up +10
            #    (3) Ages divisible by 5, greater than age at last follow up +5, and less than 80 years,
            #    to make assessments for MRI screening easier
            calc_ages = []
            alf = tage
            while alf <= settings.MAX_AGE_FOR_RISK_CALCS:
                alf += 1
                if(alf <= (tage + 5) or alf % 5 == 0 or alf == (tage + 10)):
                    calc_ages.append(alf)
        elif isinstance(calc_ages, int):
            calc_ages = [calc_ages]

        if len(calc_ages) == 0:
            calc_ages.append(0)
        if calc_ages[0] != 0:
            calc_ages.insert(0, 0)
        print("3", file=f)
        print(pedigree_file_name, file=f)
        for i, age in enumerate(calc_ages):
            print("9", file=f)
            print((age-tage if age != 0 else 0), file=f)

            print("22", file=f)
            if i < len(calc_ages)-1:
                print("yes", file=f)
            else:
                print("no", file=f)
        f.close()
        return filepath

    def write_boadicea_file(self, bwa_file=None):
        """
        Write BOADICEA pedigree file.
        """
        if bwa_file is None:
            try:
                bwa_file = open("/tmp/test.bwa.txt", "a")
            except Exception:
                return

        for person in self.people:
            print(("%-" + str(settings.MAX_LENGTH_PEDIGREE_NUMBER_STR) + "s\t") % person.famid, file=bwa_file, end="")
            print(("%-8s\t%s\t") % (person.name, person.target), file=bwa_file, end="")
            print(("%-" + str(settings.MAX_FAMILY_ID_STR_LENGTH) + "s\t" +
                   "%-" + str(settings.MAX_FAMILY_ID_STR_LENGTH) + "s\t" +
                   "%-" + str(settings.MAX_FAMILY_ID_STR_LENGTH) + "s\t") %
                  (person.pid, person.fathid, person.mothid), file=bwa_file, end="")

            print('\t'.join([person.sex(), person.mztwin, person.dead]) + '\t', file=bwa_file, end="")
            print("%-3s\t" % person.age, file=bwa_file, end="")
            print("%-4s\t" % person.yob, file=bwa_file, end="")

            d = person.cancers.diagnoses
            print("%-3s\t%-3s\t%-3s\t%-3s\t%-3s\t" %
                  (d.bc1.age if d.bc1.age != "-1" else "0",
                   d.bc2.age if d.bc2.age != "-1" else "0",
                   d.oc.age if d.oc.age != "-1" else "0",
                   d.prc.age if d.prc.age != "-1" else "0",
                   d.pac.age if d.pac.age != "-1" else "0"),
                  file=bwa_file, end="")

            record = [person.ashkn]
            # genetic tests
            gt = person.gtests
            # NOTE: order is different to settings.BC_MODEL['GENES'] so use column header
            idx = type(self).get_column_idx('Ashkn') + 1
            columns = self.get_columns()
            ngene_test_columns = 5*2
            for i in range(idx, idx+ngene_test_columns, 2):
                gene = columns[i][:-1].lower()
                gtest = getattr(gt, gene)
                record.extend([gtest.test_type, gtest.result])

            # pathology
            record.extend([p.result for p in person.pathology])

            print('\t'.join(record), file=bwa_file)
        bwa_file.flush()
        return bwa_file

    def write_boadicea_file_header(self, bwa_file=None):
        """
        Write header for BOADICEA pedigree file.
        """
        if bwa_file is None:
            try:
                bwa_file = open("/tmp/test.bwa.txt", "a")
            except Exception:
                return
        print("BOADICEA import pedigree file format 4.0", file=bwa_file)
        print("\t".join(self.get_columns()), file=bwa_file)
        bwa_file.flush()
        return bwa_file

    @classmethod
    def get_column_idx(cls, name):
        """
        Get the BOADICEA file column index from the column name
        """
        for idx, val in enumerate(cls.COLUMNS):
            if val == name or val.lower == name.lower():
                return idx
        return -1

    def get_columns(self):
        return type(self).COLUMNS


class BwaPedigree(Pedigree):
    """
    BOADICEA pedigree
    """
    COLUMNS = ["FamID", "Name", "Target", "IndivID", "FathID", "MothID", "Sex", "MZtwin", "Dead", "Age", "Yob",
               "1stBrCa", "2ndBrCa", "OvCa", "ProCa", "PanCa", "Ashkn",
               "BRCA1t", "BRCA1r", "BRCA2t", "BRCA2r", "PALB2t", "PALB2r", "ATMt", "ATMr", "CHEK2t", "CHEK2r",
               "ER", "PR", "HER2", "CK14", "CK56"]


class CanRiskPedigree(Pedigree):
    """
    CanRisk data file, contains the pedigree and optionally risk factors and PRS.
    """
    COLUMNS1 = ["FamID", "Name", "Target", "IndivID", "FathID", "MothID", "Sex", "MZtwin", "Dead", "Age", "Yob",
                "BC1", "BC2", "OC", "PRO", "PAN", "Ashkn",
                "BRCA1", "BRCA2", "PALB2", "ATM", "CHEK2", "RAD51D", "RAD51C", "BRIP1", "ER:PR:HER2:CK14:CK56"]

    COLUMNS2 = ["FamID", "Name", "Target", "IndivID", "FathID", "MothID", "Sex", "MZtwin", "Dead", "Age", "Yob",
                "BC1", "BC2", "OC", "PRO", "PAN", "Ashkn",
                "BRCA1", "BRCA2", "PALB2", "ATM", "CHEK2", "BARD1", "RAD51D", "RAD51C", "BRIP1", "ER:PR:HER2:CK14:CK56"]

    def get_prs(self, mname):
        ''' Get the PRS for the given cancer model.  '''
        if mname == 'BC' and hasattr(self, 'bc_prs'):
            return self.bc_prs
        elif mname == 'OC' and hasattr(self, 'oc_prs'):
            return self.oc_prs
        return None

    def get_rfcode(self, mname):
        ''' Get the risk factor code for the given cancer model.  '''
        if mname == 'BC' and hasattr(self, 'bc_risk_factor_code'):
            return self.bc_risk_factor_code
        elif mname == 'OC' and hasattr(self, 'oc_risk_factor_code'):
            return self.oc_risk_factor_code
        return 0

    @classmethod
    def get_column_idx(cls, name, file_type="canrisk2"):
        """
        Get the BOADICEA file column index from the column name
        """
        if file_type == "canrisk1":
            cols = cls.COLUMNS1
        else:
            cols = cls.COLUMNS2
        for idx, val in enumerate(cols):
            if val == name or val.lower == name.lower():
                return idx
        return -1
