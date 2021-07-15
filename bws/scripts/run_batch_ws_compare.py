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
from bws.pedigree import PedigreeFile

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


irates = "BOADICEA-Model-V6/Data/incidences_"+args.cancer_rates+".nml"
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


def run_batch(cwd, csvfile, ofile, irates, ashkn=False, model='BC'):
    ''' Run batch processing script. '''
    setting = FORTRAN+"settings_"+model+("_AJ" if ashkn else "")+".ini"
    cmd = [FORTRAN+"batch_run.sh",
           "-r", ofile,
           "-i", irates,
           "-s", setting,
           "-l", os.path.join(cwd, model+"runlog.log")]
    if model == 'OC':
        cmd.append('-o')
    cmd.append(csvfile)

    process = Popen(cmd, cwd=FORTRAN, stdout=PIPE, stderr=PIPE)
    (outs, errs) = process.communicate()
    _exit_code = process.wait()
    return outs, errs


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
        with open(bwa, 'r') as f:
            pedigree_data = f.read()
        pedigree = PedigreeFile(pedigree_data).pedigrees[0]

        BC_BATCH_RESULT = os.path.join(cwd, "batch_boadicea_result.out")
        outs, errs = run_batch(cwd, csvfile, BC_BATCH_RESULT, irates, ashkn=pedigree.is_ashkn())
        OC_BATCH_RESULT = os.path.join(cwd, "batch_ovarian_result.out")
        outs, errs = run_batch(cwd, csvfile, OC_BATCH_RESULT, irates, ashkn=pedigree.is_ashkn(), model='OC')

        # compare webservice.tab with batch_result.out
        f = open(args.tab, "r")
        model = 'BC'
        for line in f:
            if 'boadicea model' in line:
                model = 'BC'
            elif 'ovarian model' in line:
                model = 'OC'
            crisks = re.match("^(.*\t.+\t80\t(\d*\.\d+)).*", line)
            if crisks:
                if model == 'BC':
                    bc_80_ws = crisks.group(2)
                elif model == 'OC':
                    oc_80_ws = crisks.group(2)
        f.close()

        def get_80(fname):
            if not os.path.isfile(fname):
                print(outs)
                print(errs)

            f = open(fname, "r")
            for line in f:
                crisks = line.split(',')
                if len(crisks) == 4 and not line.startswith('file'):
                    c_80_batch = crisks[3]
            f.close()
            return c_80_batch

        bc_80_batch = get_80(BC_BATCH_RESULT)
        oc_80_batch = get_80(OC_BATCH_RESULT)

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
