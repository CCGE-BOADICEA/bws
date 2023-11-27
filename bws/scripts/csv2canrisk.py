'''

EXPERIMENTAL SCRIPT - USE WITH CAUTION

Utility script for creating CanRisk pedigree file from a CSV with
columns defined in CSV_COLUMNS

Usage:
python3 csv2canrisk.py /home/tim/ped.csv

© 2023 University of Cambridge
SPDX-FileCopyrightText: 2023 University of Cambridge
SPDX-License-Identifier: GPL-3.0-or-later
'''
import sys
import os
import csv
import argparse


RISK_FACTORS = [
    'Menarche', 'Parity', 'First_live_birth', 'OC_Use', 'OC_Duration', 'MHT_Use',
    'BMI', 'Alcohol', 'Menopause', 'BIRADS', 'Height', 'Tubal_Ligation',
    'Endometriosis', 'BC_PRS_alpha', 'BC_PRS_z', 'OC_PRS_alpha', 'OC_PRS_z']


CSV_COLUMNS = [
    "FamID", "IndivID", "FathID", "MothID", "Proband", "Sex", "MZtwin", "Age", "Yob",
    "BrCa_1st", "BrCa_2nd", "OvCa", "ProCa", "PanCa", "Censoring_Age",
    "BRCA1t", "BRCA1r", "BRCA2t", "BRCA2r", "PALB2t", "PALB2r", "CHEK2t", "CHEK2r",
    "ATMt", "ATMr", "RAD51Dt", "RAD51Dr", "RAD51Ct", "RAD51Cr",
    "BRIP1t", "BRIP1r", "BARD1t", "BARD1r", "HOXB13t", "HOXB13r",
    "ER", "PR", "HER2", "CK14", "CK56",
    "Menarche", "Parity", "First_live_birth", "OC_Use", "OC_Duration", "MHT_Use",
    "BMI", "Alcohol", "Menopause", "Mamm_density", "Height", "Tubal_Ligation", "Endometriosis",
    "BC_PRS_alpha", "BC_PRS_z", "OC_PRS_alpha", "OC_PRS_z", "PC_PRS_alpha", "PC_PRS_z"]

CANRISK_HDR = [
    "FamID", "Name", "Target", "IndivID", "FathID", "MothID", "Sex", "MZtwin", "Dead", "Age",
    "Yob", "BC1", "BC2", "OC", "PRO", "PAN", "Ashkn", "BRCA1", "BRCA2", "PALB2", "ATM", "CHEK2",
    "BARD1", "RAD51D", "RAD51C", "BRIP1", "ER:PR:HER2:CK14:CK56"]

GENES = ["BRCA1", "BRCA2", "PALB2", "ATM", "CHEK2", "BARD1", "RAD51D", "RAD51C", "BRIP1"]



def zero4NA(v):
    return v if v != "NA" and v != "" else "0"

def convert2canrisk(csvfilename):
    '''
    Convert CSV file to CanRisk file
    @param csvfilename - CSV file
    '''

    print("##CanRisk 2.0")
    
    peeps = []
    with open(csvfilename, newline='') as csvfile:
        reader = csv.DictReader(csvfile, delimiter=',')
        for row in reader:
            if row['Proband'] == "1":
                for rf in RISK_FACTORS:
                    if rf in row and row[rf] != "NA" and row[rf] != "":
                        if rf == "OC_use":
                            if "OC_Duration" in row:
                                oc_duration = rf["OC_Duration"]
                                if oc_duration != "0":
                                    print("##"+rf+"="+row[rf])
                                else:
                                    print("##"+rf+"="+row[rf]+":"+oc_duration)
                        elif rf == "OC_Duration":
                            pass
                        elif rf == "Tubal_Ligation":
                            print("##TL="+row[rf])
                        elif rf == "Endometriosis":
                            print("##Endo="+row[rf])
                        elif "_PRS_alpha" in rf:
                            print("##PRS_"+rf[:2]+"=alpha="+row[rf]+",zscore="+row[rf[:2]+"_PRS_z"])
                        elif "_PRS_" in rf:
                            pass
                        else:
                            print("##"+rf+"="+row[rf])
                
            peep = ""
            for h in CANRISK_HDR:
                if h == "Target":
                    peep+=zero4NA(row["Proband"])+"\t"
                elif h == "BC1":
                    peep+=zero4NA(row["BrCa_1st"])+"\t"
                elif h == "BC2":
                    peep+=zero4NA(row["BrCa_2nd"])+"\t"
                elif h == "OC":
                    peep+=zero4NA(row["OvCa"])+"\t"
                elif h == "PRO":
                    peep+=zero4NA(row["ProCa"])+"\t"
                elif h == "PAN":
                    peep+=zero4NA(row["PanCa"])+"\t"
                elif h == "FathID" or h == "MothID":
                    peep+=zero4NA(row[h])+"\t"
                elif h in GENES:
                    if row[h+"t"] == "NA" or row[h+"t"] == "":
                        peep+="0:0\t"
                    else:
                        peep+=row[h+"t"]+":"+row[h+"r"]+"\t"
                elif h == "ER:PR:HER2:CK14:CK56":
                    ps = ""
                    for p in h.split(":"):
                        ps+=(row[p]+":" if row[p] != "NA" and row[p] != "" else "0:")
                    peep+=ps[:-1]
                elif h not in row:
                    print("WARNING ::: "+h)
                else:
                    peep+=zero4NA(row[h])+"\t"
            peeps.append(peep.strip())
    
    print("##"+'\t'.join(CANRISK_HDR))
    for p in peeps:
        print(p)


# command line parser
if __name__ == "__main__":
    parser = argparse.ArgumentParser('EXPERIMENTAL SCRIPT TO CONVERT CSV TO CANRISK FILES- USE WITH CAUTION')
    parser.add_argument('csv',
                        type=argparse.FileType('r'),
                        help="csv input file")

    args = parser.parse_args()
    vargs = vars(args)
    csvfilename = vargs['csv'].name

    if not os.path.isfile(csvfilename):
        print(csvfilename + " is not a file.")
        sys.exit(1)
    convert2canrisk(csvfilename)
