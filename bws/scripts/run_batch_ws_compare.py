#!/usr/bin/env python3
#
# Compare PRS
# bws/scripts/run_batch_ws_compare.py -u canrisk -p bws/tests/data/canrisk_prs.txt
#

import argparse
from boadicea.scripts.boadicea_2_csv import convert2csv
from bws.scripts.run_webservice import get_auth_token, runws
from subprocess import Popen, PIPE
import tempfile
import os
import re
import shutil

# batch fortran home
FORTRAN = "/home/MINTS/tjc29/boadicea_classic/github/BOADICEA-Model-v5.0.0/"

#
# define optional command line arguments
parser = argparse.ArgumentParser('run a risk prediction via the web-service')
parser.add_argument('-c', '--can', default='boadicea', choices=['boadicea', 'ovarian', 'both'],
                    help='Cancer risk models')

parser.add_argument('--url', default='https://canrisk.org/', help='Web-services URL')
parser.add_argument('-u', '--user', help='Username')
parser.add_argument('-p', '--ped', help='CanRisk (or BOADICEA v4) pedigree file or directory of pedigree file(s)')
parser.add_argument('-t', '--tab', help='Tab delimeted output file name')
parser.add_argument('--token', help='authentication token')

args = parser.parse_args()
args.mut_freq = 'UK'
args.cancer_rates = 'UK'
url = args.url
token = get_auth_token(args, url)
bwa = input("Pedigree (BOADICEA v4/CanRisk file) or path to directory of pedigrees: ") if args.ped is None else args.ped

try:
    cwd = tempfile.mkdtemp(prefix="canrisk_batch_")

    # run webservice
    args.tab = os.path.join(cwd, 'webservice.tab')
    runws(args, {"user_id": "end_user_id"}, bwa, ['boadicea'], token, url)

    # get risk factor names and values plus PRS
    rfsnames = []
    rfs = {}
    f = open(bwa, "r")
    for line in f:
        if line.startswith("##") and "##CanRisk" not in line:
            if "PRS_BC" in line:
                zscore = re.match("##PRS.*(zscore=([-]?\d*\.\d+)).*", line)
                if zscore is not None:
                    rfsnames.append(['PRS_BC', 'PRS'])
                    rfs['PRS_BC'] = zscore[2]
            else:
                line = line.replace("##", "").strip().split("=")
                rfsnames.append([line[0], line[0].capitalize()])
                rfs[line[0]] = line[1]
    f.close()

    csvfile = os.path.join(cwd, "ped.csv")
    convert2csv(bwa, csvfile, rfsnames, rfs)

    # run batch script
    BATCH_RESULT = os.path.join(cwd, "batch_result.out")

    process = Popen(
        [FORTRAN+"batch_BOADICEA.sh",
         "-r", BATCH_RESULT,
         "-s", FORTRAN+"BOADICEA_settings.ini",
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
        risks = re.match("^(.*\t\d+\t80\t(\d*\.\d+)).*", line)
        if risks:
            bc_80_ws = risks[2]
    f.close()

    f = open(BATCH_RESULT, "r")
    for line in f:
        risks = re.match("^"+csvfile+",[^,]*(,\d+\.\d+){3},(\d+\.\d+).*", line)
        if risks:
            bc_80_batch = risks[2]
    f.close()

    if bc_80_ws and bc_80_batch and bc_80_ws == bc_80_batch:
        print("EXACT MATCH "+bc_80_ws+" "+bc_80_batch)
    else:
        print("NOT A MATCH!!!!! "+bc_80_ws+" "+bc_80_batch)
finally:
    # print(cwd)
    shutil.rmtree(cwd)
