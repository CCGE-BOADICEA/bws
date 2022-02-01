'''
Utility script for creating CSV files to be used with the batch script. As
input it takes a CanRisk pedigree file.

Dependencies:

pip install django==3.2
pip install djangorestframework==3.13.1

remove vcf2prs dependency in bws/settings.py as follows-
1) remove 'import vcf2prs' line
2) change the function get_alpha() to be:

def get_alpha(ref_file):
    return None

Usage:
export DJANGO_SETTINGS_MODULE=bws.settings
python3 -m bws.scripts.canrisk_2_csv pedigree_file.canrisk output.csv -c 1 5 80

'''
import sys
import os
import re
import argparse
from bws.pedigree import PedigreeFile
from bws.cancer import Genes, Cancers
from bws.exceptions import PedigreeFileError
# from django.conf import settings


REGEX_CANRISK1_PEDIGREE_FILE_HEADER = \
    re.compile("^(##CanRisk\\s1(.0)*)$")
REGEX_CANRISK2_PEDIGREE_FILE_HEADER = \
    re.compile("^(##CanRisk\\s2(.0)*)$")

RISK_FACTORS = ['Menarche', 'Parity', 'First_live_birth', 'OC_Use', 'OC_Duration', 'MHT_Use',
                'BMI', 'Alcohol', 'Menopause', 'BIRADS', 'Height', 'Tubal_Ligation',
                'Endometriosis', 'BC_PRS_alpha', 'BC_PRS_z', 'OC_PRS_alpha', 'OC_PRS_z']


def get_prs(line, cancer):
    '''
    Add PRS arrays for csv batch input file parameters.
    @param line: CanRisk PRS header e.g. PRS_BC=alpha=0.444,zscore=1.12
    @param cancer: string denoting cancer type, i.e. 'BC' or 'OC'
    '''
    zscore = re.match("##PRS.*(zscore=([-]?\d*\.\d+)).*", line)
    alpha = re.match("##PRS.*(alpha=([-]?\d*\.\d+)).*", line)
    if zscore is not None and alpha is not None:
        return [cancer+'_PRS_z', cancer+'_PRS_alpha'], [zscore.group(2), alpha.group(2)]


def get_rfs(line):
    '''  Get risk factor names and values plus PRS from CanRisk file for CSV file '''
    if "PRS_BC" in line:      # alpha=0.45,zscore=0.1234
        return get_prs(line, 'BC')
    elif "PRS_OC" in line:    # alpha=0.45,zscore=0.1234
        return get_prs(line, 'OC')
    else:
        line = line.replace("##", "").strip().split("=")
        if line[0] == 'TL':
            return 'Tubal_Ligation', line[1]
        elif line[0].upper() == 'MHT_USE':
            return 'MHT_Use', line[1]
        elif line[0] == 'height':
            return 'Height', line[1]
        elif line[0].upper() == 'OC_USE':
            if line[1] == "N" or line[1] == "C":
                return ['OC_Use', 'OC_Duration'], [line[1], '0']
            elif ":" in line[1]:
                parts = line[1].split(":")
                return ['OC_Use', 'OC_Duration'], [parts[0], parts[1]]
        elif line[0].upper() == 'BIRADS':
            return 'BIRADS', line[1]
        elif line[0] == 'endo':
            return 'Endometriosis', line[1]
        elif line[0].upper() == 'BMI':
            return 'BMI', line[1]
        else:
            if not line[0].capitalize() in RISK_FACTORS:
                raise PedigreeFileError(line[0].capitalize() + "NOT FOUND")
            return line[0].capitalize(), line[1]


def get_rf_values(pedigree_data):
    ''' Return risk factors as a dictionary of the FAMID. '''
    canrisk_headers = {}
    canrisk_header = {}
    famid = None
    for idx, line in enumerate(pedigree_data.splitlines()):
        if idx == 0:
            if REGEX_CANRISK1_PEDIGREE_FILE_HEADER.match(line):
                file_type = 'canrisk1'
            elif REGEX_CANRISK2_PEDIGREE_FILE_HEADER.match(line):
                file_type = 'canrisk2'
            else:
                raise PedigreeFileError(
                    "The first header record in the pedigree file has unexpected characters. " +
                    "The first header record must be '##CanRisk 2.0'.")
        elif "##CanRisk" in line or "##FamID" in line:
            pass
        elif line.startswith("##"):
            if '=' in line:                    # risk factor declaration line
                rf_name, rf_value = get_rfs(line)
                if isinstance(rf_name, list):
                    for i in range(len(rf_name)):
                        canrisk_header[rf_name[i]] = rf_value[i]
                else:
                    canrisk_header[rf_name] = rf_value
        else:
            record = line.split()
            this_famid = record[0].strip()
            this_famid = this_famid.replace("-", "")[:8]
            if famid is None or famid != this_famid:         # start of pedigree
                canrisk_headers[this_famid] = canrisk_header
                famid = this_famid
                canrisk_header = {}

            if file_type == 'canrisk1' and len(record) != 26:
                raise PedigreeFileError("A data record has an unexpected number of data items. " +
                                        "CanRisk format 1 pedigree files should have " +
                                        "26 data items per line.")
            elif file_type == 'canrisk2' and len(record) != 27:
                raise PedigreeFileError("A data record has an unexpected number of data items. " +
                                        "CanRisk format 2 pedigree files should have " +
                                        "27 data items per line.")
    return canrisk_headers


def convert2csv(filename, csvfilename, censoring_ages_freq=[1, 5, 10]):
    '''
    Convert pedigree file to a CSV file for batch processing
    @param filename - pedigree file name
    @param csvfilename - name of output CSV file
    @param censoring_ages - object of risk factors (e.g args.height)
    '''

    with open(filename, 'r') as f:
        pedigree_data = f.read()
    f.close()

    pf = PedigreeFile(pedigree_data)
    cheaders = get_rf_values(pedigree_data)
    genes = Genes.get_all_model_genes()

    hdr = ["FamID", "Name", "Proband", "IndivID", "FathID", "MothID", "Sex", "MZtwin", "Dead",
           "Age", "Yob", "BrCa_1st", "BrCa_2nd", "OvCa", "ProCa", "PanCa", "Ashkn"]

    for gene in genes:
        hdr.extend([gene+"t", gene+"r"])
    hdr.extend(["ER", "PR", "HER2", "CK14", "CK56", "Censoring_Age"])
    hdr += RISK_FACTORS     # add risk factors to the header
    csv_file = open(csvfilename, "w")
    for i in range(len(hdr)):
        print(hdr[i], file=csv_file, end="," if i < len(hdr)-1 else "")
    print('', file=csv_file)
    for ped in pf.pedigrees:
        trgt = ped.get_target()
        tage = int(trgt.age)
        alf = tage
        calc_ages = []
        while alf <= 79:
            alf += 1
            if(alf-tage in censoring_ages_freq):
                calc_ages.append(str(alf))
        calc_ages.append("80")

        for censoring_age in calc_ages:
            for person in ped.people:
                print(person.famid+":"+censoring_age, file=csv_file, end=",")
                print(person.name, file=csv_file, end=",")
                print(person.target if person.target == "1" else "", file=csv_file, end=",")
                print(person.pid, file=csv_file, end=",")
                print(person.fathid if person.fathid != "0" else "", file=csv_file, end=",")
                print(person.mothid if person.mothid != "0" else "", file=csv_file, end=",")
                print(person.sex(), file=csv_file, end=",")
                print(person.mztwin if person.mztwin != "0" else "", file=csv_file, end=",")
                print(person.dead, file=csv_file, end=",")
                print(person.age, file=csv_file, end=",")
                print(person.yob if person.yob != "0" else "", file=csv_file, end=",")

                cancers = person.cancers
                d = cancers.diagnoses
                age = ""
                [print((getattr(d, c).age if getattr(d, c).age != 'AU' and getattr(d, c).age != '-1' else age),
                       file=csv_file, end=",") for c in Cancers.get_cancers()]
                print(person.ashkn, file=csv_file, end=",")

                gtests = person.gtests
                for g in genes:
                    try:
                        gt = getattr(gtests, g.lower())
                        print(gt.test_type if gt.test_type != "0" else "", file=csv_file, end=",")
                        print(gt.result if gt.result != "0" else "", file=csv_file, end=",")
                    except AttributeError:
                        raise

                p = person.pathology
                print(p.er.result if p.er.result != "0" else "", file=csv_file, end=",")
                print(p.pr.result if p.pr.result != "0" else "", file=csv_file, end=",")
                print(p.her2.result if p.her2.result != "0" else "", file=csv_file, end=",")
                print(p.ck14.result if p.ck14.result != "0" else "", file=csv_file, end=",")
                print(p.ck56.result if p.ck56.result != "0" else "", file=csv_file, end=",")

                print((censoring_age if person.target == "1" else ""), file=csv_file, end="")
                this_rfs = None
                if person.famid in cheaders:
                    this_rfs = cheaders[person.famid]
                for rf in RISK_FACTORS:
                    if person.target == "1":
                        if this_rfs is not None:
                            if rf in this_rfs:
                                print(","+this_rfs[rf], file=csv_file, end="")
                            else:
                                print(",", file=csv_file, end="")
                        else:
                            print(",", file=csv_file, end="")
                    else:
                        print(",", file=csv_file, end="")

                print('', file=csv_file)
    csv_file.close()


# command line parser
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("pedigree", type=argparse.FileType('r'),
                        help="boadicea pedigree file")
    parser.add_argument("csv", type=argparse.FileType('w'),
                        help="csv output file")
    parser.add_argument("-c", nargs="+", type=int,
                        default=[1, 5, 10],
                        help="censoring ages")

    args = parser.parse_args()
    vargs = vars(args)
    filename = vargs['pedigree'].name
    csvfilename = vargs['csv'].name
    cen = list(args.c)

    if not os.path.isfile(filename):
        print(filename + " is not a file.")
        sys.exit(1)
    convert2csv(filename, csvfilename, cen)
