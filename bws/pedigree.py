""" Pedigree data """
import re

from django.conf import settings

from bws.cancer import Cancer, GeneticTest, PathologyTests, PathologyTest, Cancers,\
    BWSGeneticTests, CanRiskGeneticTests
from bws.exceptions import PedigreeFileError, PedigreeError, PersonError
from datetime import date
from random import randint
import abc
from bws.risk_factors.bc import BCRiskFactors
from bws.risk_factors.oc import OCRiskFactors


# BOADICEA header
REGEX_BWA_PEDIGREE_FILE_HEADER_ONE = \
    re.compile("^(BOADICEA\\s+import\\s+pedigree\\s+file\\s+format\\s[124](.0)*)$")
REGEX_CANRISK_PEDIGREE_FILE_HEADER_ONE = \
    re.compile("^(##CanRisk\\s[1](.0)*)$")
REGEX_ALPHANUM_HYPHENS = re.compile("^([\\w\-]+)$")
REGEX_ONLY_HYPHENS = re.compile("^([\-]+)$")
REGEX_ONLY_ZEROS = re.compile("^[0]+$")
REGEX_AGE = re.compile("^\\d{1,3}$")
REGEX_YEAR_OF_BIRTH = re.compile("^0|((17|18|19|20)[0-9][0-9])$")
REGEX_ASHKENAZI_STATUS = re.compile("^[01]$")

BLANK_LINE = re.compile(r'^\s*$')

# calculation type
CANCER_RISKS = 1
MUTATION_PROBS = 2


class CanRiskHeader():
    '''
    CanRisk File Format Header
    '''

    def __init__(self):
        self.lines = []

    def add_line(self, line):
        ''' Append header line to the list of lines. '''
        self.lines.append(line)

    def get_risk_factor_codes(self):
        ''' Get risk factor code from header lines. '''
        bc_risk_factors = BCRiskFactors()
        oc_risk_factors = OCRiskFactors()
        for line in self.lines:
            try:
                parts = line.split('=')
                rfnam = parts[0][2:].lower()    # risk factor name
                rfval = parts[1]                # risk factor value

                # lookup breast/ovarian cancer risk factors
                bc_risk_factors.add_category(rfnam, rfval)
                oc_risk_factors.add_category(rfnam, rfval)
            except Exception as e:
                raise PedigreeFileError(e)
        return BCRiskFactors.encode(bc_risk_factors.cats), OCRiskFactors.encode(oc_risk_factors.cats)


class PedigreeFile(object):
    """
    CanRisk and BOADICEA import pedigree file.
    """

    def __init__(self, pedigree_data):
        self.pedigree_data = pedigree_data
        pedigrees_records = [[]]
        canrisk_headers = []
        canrisk_header = CanRiskHeader()
        pid = 0
        famid = None
        file_type = None

        for idx, line in enumerate(pedigree_data.splitlines()):
            if idx == 0:
                if REGEX_CANRISK_PEDIGREE_FILE_HEADER_ONE.match(line):
                    file_type = 'canrisk'
                elif REGEX_BWA_PEDIGREE_FILE_HEADER_ONE.match(line):
                    file_type = 'bwa'
                else:
                    raise PedigreeFileError(
                        "The first header record in the pedigree file has unexpected characters. " +
                        "The first header record must be 'BOADICEA import pedigree file format 4'.")
            elif line.startswith('##'):
                if '=' in line:                     # risk factor declaration line
                    canrisk_header.add_line(line)
            elif idx == 1:
                self.column_names = line.split()
                if (((self.column_names[0] != 'FamID') or
                     (self.column_names[2] != 'Target') or
                     (self.column_names[3] != 'IndivID') or
                     (self.column_names[4] != 'FathID') or
                     (self.column_names[5] != 'MothID'))):
                    raise PedigreeFileError(
                        "Column headers in the pedigree file contains unexpected characters. " +
                        "It must include the 'FamID', 'Target', 'IndivID','FathID' and 'MothID' " +
                        "in columns 1, 3, 4, 5 and 6 respectively.")
            elif BLANK_LINE.match(line):
                continue
            else:
                record = line.split()
                if famid is None or famid != record[0]:         # start of pedigree
                    canrisk_headers.append(canrisk_header)
                    canrisk_header = CanRiskHeader()
                if famid is not None and famid != record[0]:    # start of next pedigree found
                    pedigrees_records.append([])
                    pid += 1
                famid = record[0]

                if(file_type == 'bwa' and len(record) != settings.BOADICEA_PEDIGREE_FORMAT_FOUR_DATA_FIELDS):
                    raise PedigreeFileError("A data record has an unexpected number of data items. " +
                                            "BOADICEA format 4 pedigree files should have " +
                                            str(settings.BOADICEA_PEDIGREE_FORMAT_FOUR_DATA_FIELDS) +
                                            " data items per line.")
                elif(file_type == 'canrisk' and len(record) != settings.BOADICEA_CANRISK_FORMAT_ONE_DATA_FIELDS):
                    raise PedigreeFileError("A data record has an unexpected number of data items. " +
                                            "CanRisk format 1 pedigree files should have " +
                                            str(settings.BOADICEA_CANRISK_FORMAT_ONE_DATA_FIELDS) +
                                            " data items per line.")
                pedigrees_records[pid].append(line)

        self.pedigrees = []
        for i in range(pid+1):
            if file_type == 'bwa':
                self.pedigrees.append(BwaPedigree(pedigree_records=pedigrees_records[i], file_type=file_type))
            elif file_type == 'canrisk':
                try:
                    bc_rfc, oc_rfc = canrisk_headers[i].get_risk_factor_codes()
                except Exception:
                    bc_rfc = 0
                    oc_rfc = 0
                self.pedigrees.append(CanRiskPedigree(pedigree_records=pedigrees_records[i], file_type=file_type,
                                                      bc_risk_factor_code=bc_rfc, oc_risk_factor_code=oc_rfc))

    @classmethod
    def validate(cls, pedigrees):
        if isinstance(pedigrees, Pedigree):
            pedigrees = [pedigrees]
        warnings = []
        for pedigree in pedigrees:
            people = pedigree.people
            pedigree.validate()             # Validate pedigree input data

            for p in people:
                if not p.is_complete():
                    warnings.append("year of birth and age at last follow up must be specified in order for " +
                                    p.pid + " to be included in a calculation")
                p.validate(pedigree)                        # Validate person data
                type(p.cancers).validate(p)                 # Validate cancer diagnoses
                warnings.extend(PathologyTest.validate(p))  # Validate pathology status
                GeneticTest.validate(p)                     # Validate genetic tests

        return warnings


class Pedigree(metaclass=abc.ABCMeta):
    """
    A pedigree object.
    """

    def __init__(self, pedigree_records=None, people=None, file_type=None,
                 bc_risk_factor_code=None, oc_risk_factor_code=None):
        """
        @keyword pedigree_records: the pedigree records section of the BOADICEA import pedigree file.
        @keyword people: members of the pedigree.
        @keyword file_type: file type is 'bwa' or 'canrisk'.
        """
        self.people = []
        if pedigree_records is not None:
            self.famid = pedigree_records[0].split()[0]
            ids = []
            for record in pedigree_records:
                p = Person.factory(record, file_type=file_type)
                if p.target != '0' and p.target != '1':
                    raise PedigreeError("A value in the Target data column has been set to '" + p.target +
                                        "'. Target column parameters must be set to '0' or '1'.")
                if p.is_target():
                    self.target = p

                if p.pid in ids:
                    raise PedigreeError("Individual ID '" + p.pid + "' appears more than once in the pedigree file.")
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
                                "index individuals. Only one target can be specified.")
        if pedigree_size > settings.MAX_PEDIGREE_SIZE or pedigree_size < settings.MIN_BASELINE_PEDIGREE_SIZE:
            raise PedigreeError("Pedigree (" + self.famid + ") has unexpected number of family members " +
                                str(pedigree_size))
        if file_type == 'canrisk':
            if bc_risk_factor_code is not None:
                self.bc_risk_factor_code = bc_risk_factor_code
            if oc_risk_factor_code is not None:
                self.oc_risk_factor_code = oc_risk_factor_code

    def validate(self):
        """ Validation check for pedigree input.
        @param p: Person to validate pedigree data.
        """
        if(len(self.famid) > settings.MAX_LENGTH_PEDIGREE_NUMBER_STR or
           not REGEX_ALPHANUM_HYPHENS.match(self.famid) or       # must be alphanumeric plus hyphen
           REGEX_ONLY_HYPHENS.match(self.famid) or               # but not just hyphens
           REGEX_ONLY_ZEROS.match(self.famid)):                  # and not just zeros
            raise PedigreeError(
                "Family ID (1st data column) has been set to '" + self.famid +
                "'. Family IDs must be specified with between 1 and "+str(settings.MAX_LENGTH_PEDIGREE_NUMBER_STR) +
                " non-zero number or alphanumeric characters.")

        unconnected = self.unconnected()
        if len(unconnected) > 0:
            raise PedigreeError("Pedigree (" + self.famid + ") family members are not physically " +
                                "connected to the target: " + str(unconnected))

        # Check that the index's parameters are valid
        target = self.get_target()
        if target.yob == '0':
            raise PedigreeError("The target's year of birth has been set to '" + target.yob +
                                "'. This person must be assigned a valid year of birth.")
        if target.age == '0':
            raise PedigreeError("The target's age has been set to '" + target.age +
                                "'. This person must be assigned an age.")

        # Check that carrier probabilities / cancer risks can be computed
        carrier_probs = self.is_carrier_probs_viable(target=target)
        cancer_risks = self.is_risks_calc_viable(target=target)
        if(not carrier_probs and not cancer_risks):
            raise PedigreeError(
                "BOADICEA cannot compute mutation carrier probabilities because the target '" + target.pid +
                "' has a positive genetic test. Also BOADICEA cannot compute breast and ovarian cancer "
                "risks because the target is: (1) over " + str(settings.MAX_AGE_FOR_RISK_CALCS) +
                " years old or (2) male, or (3) an affected female who has developed contralateral "
                "breast cancer, ovarian cancer or pancreatic cancer.")

        #
        # Check monozygotic (MZ) twin data
        twin_store = self.get_twins()

        # Check that MZ siblings are only specified as twins, no identical triplets etc
        for t in twin_store:
            twins = twin_store[t]
            if len(twins) != 2:
                raise PedigreeError(
                    "MZ twin identifier '" + str(twins[0].pid) + "' does not appear twice in the pedigree file. "
                    "Only MZ twins are permitted in the pedigree, MZ triplets or quads are not allowed.")

            # Check MZ twin characters are valid
            if len(t) != 1 or t not in settings.UNIQUE_TWIN_IDS:
                raise PedigreeError("Invalid MZ twin character '" + t + "'. MZ twins must be identified using one " +
                                    "of the following ASCII characters: " + str(settings.UNIQUE_TWIN_IDS) + ".")

            # Check that monozygotic (MZ) twin data are consistent
            if(twins[0].mothid != twins[1].mothid or
               twins[0].fathid != twins[1].fathid):
                raise PedigreeError("Monozygotic (MZ) twins identified with the character '" + t + "' have different "
                                    "parents. MZ twins must have the same parents.")
            if(twins[0].yob != twins[1].yob):
                raise PedigreeError("Monozygotic (MZ) twins identified with the character '" + t + "' have different "
                                    "years of birth. MZ twins must have the same year of birth.")

            # Check that living MZ twins have the same age at last follow up
            if(twins[0].dead == '0' and twins[1].dead == '0' and twins[0].age != twins[1].age):
                raise PedigreeError("Monozygotic (MZ) twins identified with the character '" + t + "' have different "
                                    "ages. If both MZ twins are alive, they must have the same age at last follow up.")

            if twins[0].sex() != twins[1].sex():
                raise PedigreeError("Monozygotic (MZ) twins identified with the character '" + t + "' have a different "
                                    "sex. MZ twins must have the same sex.")

            # Check that the MZ twins have the same genetic status
            if not GeneticTest.compareTestResults(twins[0], twins[1]):
                raise PedigreeError("Monozygotic (MZ) twins have both had a genetic test, but the genetic test results "
                                    "for these individuals are different. Under these circumstances, the genetic test "
                                    "results must be the same.")

        # Check to ensure that the maximum number of MZ twin pairs per pedigree has not been exceeded
        if len(twin_store.keys()) > settings.MAX_NUMBER_MZ_TWIN_PAIRS:
            raise PedigreeError("Maximum number of MZ twin pairs has been exceeded. Input pedigrees must have a "
                                "maximum of " + str(settings.MAX_NUMBER_MZ_TWIN_PAIRS) + " MZ twin pairs.")

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
             (d.bc2.age != "0") or (d.pac.age != "0") or (d.oc.age != "0") or
             (int(target.age) > settings.MAX_AGE_FOR_RISK_CALCS) or
             (int(target.yob) < settings.MIN_YEAR_OF_BIRTH))):
            return False
        return True

    def is_carrier_probs_viable(self,  target=None):
        """
        Return true if the target does not have a positive genetic test carrier probs
        cannot be calculated.
        @return: true if carrier probability calculation is viable
        """
        if target is None:
            target = self.get_target()
        gtests = target.gtests
        for t in gtests:
            if t.result == 'P':
                return False
        return True

    def write_pedigree_file(self, file_type, risk_factor_code='0', prs=None, filepath="/tmp/test.ped",
                            model_settings=settings.BC_MODEL):
        """
        Write input pedigree file for fortran.
        """
        f = open(filepath, "w")
        print("(I3,X,A8)", file=f)
        if file_type == MUTATION_PROBS:
            pcount = (len(model_settings['GENES'])+1)
        else:
            if model_settings['NAME'] == "BC":
                pcount = 3
            elif model_settings['NAME'] == "OC":
                pcount = 2

        risk_factor_str_len = str(len(risk_factor_code))

        print("(3(A7,X),2(A1,X),2(A3,X)," + str(len(model_settings['CANCERS'])+1) + "(A3,X),5(A1,X),A4," +
              ("6(X,A1),X,A"+risk_factor_str_len+",2(X,A6))"
               if (model_settings['NAME'] == 'BC') else "1(X,A1),3(X,A8))"),
              file=f)

        for gt in range(pcount):
            print("%-3d %-8s" % (len(self.people), self.people[0].famid), file=f)

            for p in self.people:
                # IndivID FathID MothID Sex MZ Genotype, Polygene 1BC 2BC OC
                genotype = gt if (p.target != "0" and file_type == MUTATION_PROBS) else ''
                proband_status = (gt+1) if (p.target != "0" and file_type == CANCER_RISKS) else p.target
                print("%-7s %-7s %-7s %-1s %-1s %3s %-3s " %
                      (p.pid,
                       p.fathid if p.fathid != "0" else '',
                       p.mothid if p.mothid != "0" else '', p.sex(),
                       p.mztwin if p.mztwin != "0" else '',
                       genotype, '   '), file=f, end="")

                print(p.cancers.write(model_settings['CANCERS']), file=f, end="")
                print("%3s " % p.age, file=f, end="")

                # Gene Tests
                gtests = p.gtests
                for g in model_settings['GENES']:
                    print("%1s " % getattr(gtests, g.lower()).get_genetic_test_data(), file=f, end="")

                print("%4s " % (p.yob if p.yob != 0 else settings.MENDEL_NULL_YEAR_OF_BIRTH), file=f, end="")

                if model_settings['NAME'] == 'BC':
                    print(PathologyTest.write(p.pathology), file=f, end="")

                # ProbandStatus RiskFactor PolygStanDev PolygLoad
                if model_settings['NAME'] == 'BC':
                    fmt = "%1s %"+risk_factor_str_len+"s %6.5f %6.5f"
                else:
                    fmt = "%1s %8s %8.4f %8.4f"
                print(fmt % (proband_status, risk_factor_code,
                             prs.alpha if prs is not None and prs.alpha else 0,
                             prs.beta if prs is not None and prs.beta else 0,), file=f)
        f.close()
        return filepath

    def write_batch_file(self, batch_type, pedigree_file_name, filepath="/tmp/test.bat",
                         mutation_freq=settings.BC_MODEL['MUTATION_FREQUENCIES']['UK'],
                         sensitivity=settings.BC_MODEL['GENETIC_TEST_SENSITIVITY'],
                         calc_ages=None):
        """
        Write fortran input batch file.
        @param batch_type: compute MUTATION_PROBS or CANCER_RISKS
        @param pedigree_file_name: path to fortran pedigree file
        @param filepath: path to write the batch file to
        @param mutation_freq: mutation frequencies
        @param sensitivity: genetic test sensitivity
        @param calc_ages: list of ages to calculate a cancer risk at
        """
        f = open(filepath, "w")

        if (batch_type != MUTATION_PROBS) and (batch_type != CANCER_RISKS):
            raise PedigreeFileError(
                    "Invalid batch file type.")

        if 'PALB2' in mutation_freq:
            model_settings = settings.BC_MODEL
        elif 'RAD51D' in mutation_freq:
            model_settings = settings.OC_MODEL

        if batch_type == MUTATION_PROBS:
            print("3", file=f)
            print(pedigree_file_name, file=f)
            print("8", file=f)
            for gene in model_settings['GENES']:
                print(mutation_freq[gene], file=f)
            for gene in model_settings['GENES']:
                print(sensitivity[gene], file=f)

            print("21", file=f)
            print("no", file=f)
        elif batch_type == CANCER_RISKS:
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

            print("3", file=f)
            print(pedigree_file_name, file=f)
            for i, age in enumerate(calc_ages):
                print("8", file=f)
                print((age-tage), file=f)
                for gene in model_settings['GENES']:
                    print(mutation_freq[gene], file=f)
                for gene in model_settings['GENES']:
                    print(sensitivity[gene], file=f)

                print("21", file=f)
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
                  (d.bc1.age, d.bc2.age, d.oc.age, d.prc.age, d.pac.age),
                  file=bwa_file, end="")

            record = [person.ashkn]
            # genetic tests
            gt = person.gtests
            # NOTE: order is different to settings.BC_MODEL['GENES'] so use column header
            idx = type(self).get_column_idx('Ashkn') + 1
            columns = self.get_columns()
            ngene_test_columns = len(settings.BC_MODEL['GENES'])*2
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
    CanRisk pedigree
    """
    COLUMNS = ["FamID", "Name", "Target", "IndivID", "FathID", "MothID", "Sex", "MZtwin", "Dead", "Age", "Yob",
               "BC1", "BC2", "OC", "PRO", "PAN", "Ashkn",
               "BRCA1", "BRCA2", "PALB2", "ATM", "CHEK2", "RAD51D", "RAD51C", "BRIP1", "ER:PR:HER2:CK14:CK56"]


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
           not REGEX_ALPHANUM_HYPHENS.match(self.name)):
            raise PersonError("A name '"+self.name+"' is unspecified or is not an alphanumeric string.")

        if(len(self.pid) < settings.MIN_FAMILY_ID_STR_LENGTH or
           len(self.pid) > settings.MAX_FAMILY_ID_STR_LENGTH or
           REGEX_ONLY_ZEROS.match(self.pid) or
           not REGEX_ALPHANUM_HYPHENS.match(self.pid)):
            raise PersonError("An individual identifier (IndivID column) was specified as '" + self.pid +
                              ". Individual identifiers must be alphanumeric strings with a maximum of " +
                              str(settings.MAX_FAMILY_ID_STR_LENGTH)+"characters.")

        if(len(self.fathid) < settings.MIN_FAMILY_ID_STR_LENGTH or
           len(self.fathid) > settings.MAX_FAMILY_ID_STR_LENGTH or
           not REGEX_ALPHANUM_HYPHENS.match(self.fathid)):
            raise PersonError("Father identifier ('" + self.fathid + "', FathID column) has unexpected characters. "
                              "It must be alphanumeric strings with a maximum of " +
                              str(settings.MAX_FAMILY_ID_STR_LENGTH) + " characters")

        if(len(self.mothid) < settings.MIN_FAMILY_ID_STR_LENGTH or
           len(self.mothid) > settings.MAX_FAMILY_ID_STR_LENGTH or
           not REGEX_ALPHANUM_HYPHENS.match(self.mothid)):
            raise PersonError("Mother identifier ('" + self.mothid + "', MothID column) has unexpected characters. "
                              "It must be alphanumeric strings with a maximum of " +
                              str(settings.MAX_FAMILY_ID_STR_LENGTH) + " characters")

        if(self.fathid == '0' and self.mothid != '0') or (self.fathid != '0' and self.mothid == '0'):
            raise PersonError("Family member '"+self.name+"' has only one parent specified. All family members must "
                              "have no parents specified (i.e. they must be founders) or both parents specified.")

        # check for missing parents
        if self.mothid != '0' and pedigree.get_person(self.mothid) is None:
            raise PersonError("The mother '"+self.mothid+"' of family member '" + self.pid +
                              "' is missing from the pedigree.")
        elif self.fathid != '0' and pedigree.get_person(self.fathid) is None:
            raise PersonError("The father '"+self.fathid+"' of family member '" + self.pid +
                              "' is missing from the pedigree.")

        # check all fathers are male
        if self.fathid != '0' and pedigree.get_person(self.fathid).sex() != 'M':
            raise PersonError("The father of family member '" + self.pid + "' is not specified as male. " +
                              "All fathers in the pedigree must have sex specified as 'M'.")
        # check all mothers are female
        if self.mothid != '0' and pedigree.get_person(self.mothid).sex() != 'F':
            raise PersonError("The mother of family member '" + self.pid + "' is not specified as female. " +
                              "All mothers in the pedigree must have sex specified as 'F'.")
        # check if the dead attribute has been correctly set
        if self.dead != '0' and self.dead != '1':
            raise PersonError("The family member '" + self.pid + "' has an invalid vital status " +
                              "(alive must be specified as '0', and dead specified as '1')")
        # check that age of last follow up set to either 0 (unknown) or in range 1-110
        if not REGEX_AGE.match(self.age) or int(self.age) > settings.MAX_AGE:
            raise PersonError("The age specified for family member '" + self.pid + "' has unexpected " +
                              "characters. Ages must be specified with as '0' for unknown, or in the " +
                              "range 1-" + str(settings.MAX_AGE))

        # validate year of birth
        current_year = date.today().year
        if self.yob != "0":
            if(not REGEX_YEAR_OF_BIRTH.match(self.yob) or
               int(self.yob) < settings.MIN_YEAR_OF_BIRTH or
               int(self.yob) > current_year):
                raise PersonError("The year of birth '" + self.yob + "' specified for family member '" + self.pid +
                                  "' is out of range. Years of birth must be in the range " +
                                  str(settings.MIN_YEAR_OF_BIRTH) + "-" + str(current_year))

        # Check that we have valid data values for the Ashkenazi flag
        if(not REGEX_ASHKENAZI_STATUS.match(self.ashkn)):
            raise PersonError("Family member '" + self.pid + "' has been assigned an invalid Ashkenazi "
                              "origin parameter. The Ashkenazi origin parameter must be set to '1' "
                              "for Ashkenazi origin, or '0' for not Ashkenazi origin.")

        # Check max no. of siblings not exceeded
        (siblings, siblings_same_yob) = pedigree.get_siblings(self)
        if len(siblings) > settings.MAX_NUMBER_OF_SIBS_PER_NUCLEAR_FAMILY:
            raise PersonError("Family member '" + self.pid + "' exeeded the maximum number of siblings (" +
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
        cancers = Cancers(bc1=Cancer(cols[11]), bc2=Cancer(cols[12]), oc=Cancer(cols[13]),
                          prc=Cancer(cols[14]), pac=Cancer(cols[15]))

        # use column headers to get gene test type and result
        if file_type == 'bwa':
            gtests = BWSGeneticTests.factory([GeneticTest(cols[BwaPedigree.get_column_idx(gene+'t')],
                                                          cols[BwaPedigree.get_column_idx(gene+'r')])
                                             for gene in settings.BC_MODEL['GENES']])
            pathology = PathologyTests(
                            er=PathologyTest(PathologyTest.ESTROGEN_RECEPTOR_TEST, cols[27]),
                            pr=PathologyTest(PathologyTest.PROGESTROGEN_RECEPTOR_TEST, cols[28]),
                            her2=PathologyTest(PathologyTest.HER2_TEST, cols[29]),
                            ck14=PathologyTest(PathologyTest.CK14_TEST, cols[30]),
                            ck56=PathologyTest(PathologyTest.CK56_TEST, cols[31]))
        else:
            genes = settings.BC_MODEL['GENES'] + settings.OC_MODEL['GENES'][2:]

            def get_genetic_test(cols, gene):
                idx = CanRiskPedigree.get_column_idx(gene)
                if idx < 0:
                    raise PedigreeError("Genetic test column for '" + gene + "not found.")
                gt = cols[idx].split(':')
                return GeneticTest(gt[0], gt[1])
            gtests = CanRiskGeneticTests.factory([get_genetic_test(cols, gene) for gene in genes])

            path = cols[len(CanRiskPedigree.COLUMNS)-1].split(':')
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
