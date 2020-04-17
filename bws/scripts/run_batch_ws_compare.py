#!/usr/bin/env python3
#
# Test script to compare results from BC webservice and batch script
# bws/scripts/run_batch_ws_compare.py -u USERNAME -p bws/tests/data/canrisk_prs.txt
# bws/scripts/run_batch_ws_compare.py -u USERNAME -p bws/tests/data/batch/ --url http://0.0.0.0:8000/ \
#                                     --fortran /home/xxxx/Model-Batch-Processing --cancer_rates Spain
#

import argparse
import os
from os.path import expanduser
import re
import shutil
from subprocess import Popen, PIPE
import tempfile

from boadicea.scripts.boadicea_2_csv import convert2csv
from bws.scripts.run_webservice import get_auth_token, runws
import math

#
# define optional command line arguments
parser = argparse.ArgumentParser('run a risk prediction via the web-service')

parser.add_argument('--url', default='https://canrisk.org/', help='Web-services URL')
parser.add_argument('-u', '--user', help='Username')
parser.add_argument('-p', '--ped', help='CanRisk (or BOADICEA v4) pedigree file or directory of pedigree file(s)')
parser.add_argument('-f', '--fortran', help='Path to BOADICEA model code',
                    default=os.path.join(expanduser("~"), "boadicea_classic/github/Model-Batch-Processing/"))
parser.add_argument('--cancer_rates', default='UK',
                    choices=['UK', 'Australia', 'Canada', 'USA', 'Denmark', 'Finland',
                             'Iceland', 'New-Zealand', 'Norway', 'Spain', 'Sweden'],
                    help='Cancer incidence rates (default: %(default)s)')
parser.add_argument('--token', help='authentication token')

args = parser.parse_args()
args.mut_freq = 'UK'


irates = "BOADICEA-Model/Data/incidence_rates_"+args.cancer_rates+".nml"
print('Cancer Incidence Rates: '+args.cancer_rates)

url = args.url

# batch fortran home
FORTRAN = args.fortran
if not os.path.exists(os.path.join(FORTRAN, 'batch_run.sh')):
    print("Error: check path to fortran : "+FORTRAN)
    exit(1)

token = get_auth_token(args, url)
bwa = input("Pedigree (BOADICEA v4/CanRisk file) or path to directory of pedigrees: ") if args.ped is None else args.ped
if os.path.isfile(bwa):
    bwalist = [bwa]
else:
    bwalist = [os.path.join(bwa, f) for f in os.listdir(bwa) if os.path.isfile(os.path.join(bwa, f))]


def add_prs(line, cancer, rfsnames, rfs):
    '''
    Add PRS arrays for csv batch input file parameters.
    @param line: CanRisk PRS header e.g. PRS_BC=alpha=0.444,zscore=1.12
    @param cancer: string donating cancer type, i.e. 'BC' or 'OC'
    @param rfsnames: array of risk factor names
    @param rfs: risk factor values
    '''
    zscore = re.match("##PRS.*(zscore=([-]?\d*\.\d+)).*", line)
    alpha = re.match("##PRS.*(alpha=([-]?\d*\.\d+)).*", line)
    if zscore is not None:
        rfsnames.append(['PRS_'+cancer+'_z', cancer+'_PRS_z'])
        rfs['PRS_'+cancer+'_z'] = zscore.group(2)
    if alpha is not None:
        rfsnames.append(['PRS_'+cancer+'_alpha', cancer+'_PRS_alpha'])
        rfs['PRS_'+cancer+'_alpha'] = alpha.group(2)


# loop over canrisk files and compare results from webservices with those from the batch script
exact_matches = 0
for bwa in bwalist:
    try:
        cwd = tempfile.mkdtemp(prefix="canrisk_batch_")

        # run webservice
        args.tab = os.path.join(cwd, 'webservice.tab')
        runws(args, {"user_id": "end_user_id"}, bwa, ['boadicea', 'ovarian'], token, url)

        # get risk factor names and values plus PRS
        rfsnames = []
        rfs = {}
        f = open(bwa, "r")
        for line in f:
            if line.startswith("##") and "##CanRisk" not in line and "##FamID" not in line:
                if "PRS_BC" in line:      # alpha=0.45,zscore=0.1234
                    add_prs(line, 'BC', rfsnames, rfs)
                elif "PRS_OC" in line:    # alpha=0.45,zscore=0.1234
                    add_prs(line, 'OC', rfsnames, rfs)
                else:
                    line = line.replace("##", "").strip().split("=")

                    if line[0].isupper():
                        if line[0] == 'TL':
                            name = 'Tubal_Ligation'
                        else:
                            name = line[0]
                    elif line[0] == 'mht_use':
                        name = 'MHT_use'
                    elif line[0] == 'oc_use':
                        name = 'OC_Use'
                        if line[1] == "N":
                            rfsnames.append(['OC_Duration', 'OC_Duration'])
                            rfs['OC_Duration'] = '0'
                    elif line[0] == 'birads':
                        name = 'BIRADS'
                    elif line[0] == 'endo':
                        name = 'Endometriosis'
                    else:
                        name = line[0].capitalize()

                    rfsnames.append([line[0], name])
                    rfs[line[0]] = line[1]
        f.close()

        csvfile = os.path.join(cwd, "ped.csv")
        convert2csv(bwa, csvfile, rfsnames, rfs)

        # run batch script
        BATCH_RESULT = os.path.join(cwd, "batch_result.out")

        process = Popen(
            [FORTRAN+"batch_run.sh",
             "-r", BATCH_RESULT,
             "-s", FORTRAN+"settings.ini",
             "-i", irates,
             "-l", os.path.join(cwd, "runlog.log"),
             csvfile],
            cwd=FORTRAN,
            stdout=PIPE,
            stderr=PIPE)
        (outs, errs) = process.communicate()
        exit_code = process.wait()

        # compare webservice.tab with batch_result.out
        f = open(args.tab, "r")
        for line in f:
            bcrisks = re.match("^(.*\t.+\t80\t(\d*\.\d+)).*\).*\).*\).*\)", line)
            ocrisks = re.match("^(.*\t.+\t80\t(\d*\.\d+)).*", line)
            if bcrisks:
                bc_80_ws = bcrisks.group(2)
            if ocrisks:
                oc_80_ws = ocrisks.group(2)
        f.close()

        def get_80(fname):
            f = open(fname, "r")
            for line in f:
                crisks = re.match("^"+csvfile+",[^,]*(,\d+\.\d+){3},(\d+\.\d+).*", line)
                if crisks:
                    c_80_batch = crisks.group(2)
            f.close()
            return c_80_batch

        bc_80_batch = get_80(BATCH_RESULT+"_boadicea.csv")
        oc_80_batch = get_80(BATCH_RESULT+"_ovarian.csv")

        if bc_80_ws and bc_80_batch and math.isclose(float(bc_80_ws), float(bc_80_batch)):
            print("BC EXACT MATCH ::: "+bwa+"    webservice: "+bc_80_ws+" batch: "+bc_80_batch)
        else:
            print("BC NOT A MATCH ::: "+bwa+"    webservice: "+bc_80_ws+" batch: "+bc_80_batch)
            exact_matches += 1

        if oc_80_ws and oc_80_batch and math.isclose(float(oc_80_ws), float(oc_80_batch)):
            print("OC EXACT MATCH ::: "+bwa+"    webservice: "+oc_80_ws+" batch: "+oc_80_batch)
        else:
            print("OC NOT A MATCH ::: "+bwa+"    webservice: "+oc_80_ws+" batch: "+oc_80_batch)
            exact_matches += 1
    finally:
        shutil.rmtree(cwd)

exit(exact_matches)
