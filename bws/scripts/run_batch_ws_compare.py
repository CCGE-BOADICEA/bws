#!/usr/bin/env python3
#
# Test script to compare results from BC webservice and batch script
# bws/scripts/run_batch_ws_compare.py -u USERNAME -p bws/tests/data/canrisk_prs.txt
# bws/scripts/run_batch_ws_compare.py -u USERNAME -p bws/tests/data/batch/ --url http://0.0.0.0:8000/ \
#                                     --fortran /home/xxxx/Model-Batch-Processing --cancer_rates Spain
#
# © 2022 Cambridge University
# SPDX-FileCopyrightText: 2022 Cambridge University
# SPDX-License-Identifier: GPL-3.0-or-later

from boadicea.scripts.boadicea_2_csv import convert2csv
from bws.scripts.batch import run_batch, get_censoring_ages, get_batch_results, get_rfs,\
    get_mp, compare_mp
from bws.scripts.run_webservice import get_auth_token, runws
from os.path import expanduser
import argparse
import math
import os
import re
import shutil
import tempfile


def get_ws_results(args, calc_ages):
    ''' Parse web-service tab file and return breast and ovarian cancer risks and mutation carrier probabilities. '''
    bc_ws, oc_ws, bc_mp_ws, oc_mp_ws = {}, {}, None, None
    sages = '|'.join([str(c) for c in calc_ages])
    with open(args.tab, 'r') as f:
        model = 'BC'
        for line in f:
            if 'boadicea model' in line:
                model = 'BC'
            elif 'ovarian model' in line:
                model = 'OC'
            crisks = re.match("^(.*\t.+\t("+sages+")\t(\d*\.?\d*[e\-\d*]*)).*", line)
            if crisks:
                if model == 'BC':
                    bc_ws[int(crisks.group(2))] = crisks.group(3)
                elif model == 'OC':
                    oc_ws[int(crisks.group(2))] = crisks.group(3)
            elif line.startswith('no mutation'):
                gkeys = line.strip().split('\t')
                vals = next(f).strip().split('\t')
                if model == 'BC':
                    bc_mp_ws = {gkeys[idx]: val for idx, val in enumerate(vals)}
                else:
                    oc_mp_ws = {gkeys[idx]: val for idx, val in enumerate(vals)}
        f.close()
    return bc_ws, oc_ws, bc_mp_ws, oc_mp_ws


def glob_re(path_pattern):
    ''' Get files in a path that match a pattern (e.g. /path/xyz[1-2]?). '''
    pdir = path_pattern.rsplit("/", 1)
    fnames = filter(re.compile(pdir[1]).match, os.listdir(pdir[0]))
    return [os.path.join(pdir[0], f) for f in fnames]


#
# define optional command line arguments
parser = argparse.ArgumentParser('run a risk prediction via the web-service')

parser.add_argument('--url', default='https://canrisk.org/', help='Web-services URL')
parser.add_argument('-u', '--user', help='Username')
parser.add_argument('-p', '--ped', help='CanRisk (or BOADICEA v4) pedigree file or directory of pedigree file(s)')
parser.add_argument('-f', '--fortran', help='Path to BOADICEA model code',
                    default=os.path.join(expanduser("~"), "boadicea_classic/github/Model-Batch-Processing/"))
parser.add_argument('--cancer_rates', default='UK',
                    choices=['UK', 'Australia', 'Canada', 'USA', 'Denmark', 'Estonia', 'Finland', 'France',
                             'Iceland', 'Netherlands', 'New-Zealand', 'Norway', 'Slovenia', 'Spain', 'Sweden'],
                    help='Cancer incidence rates (default: %(default)s)')
parser.add_argument('--token', help='authentication token')
parser.add_argument('--bc_rr_tolerance', default=1e-09, help='BC tolerance comparing web-service & batch risks')
parser.add_argument('--oc_rr_tolerance', default=1e-09, help='OC tolerance comparing web-service & batch risks')

parser.add_argument('--bc_probs_tolerance', default=1e-09, help='BC tolerance comparing web-service & batch probs')
parser.add_argument('--oc_probs_tolerance', default=1e-09, help='OC tolerance comparing web-service & batch probs')

args = parser.parse_args()
args.mut_freq = 'UK'

bc_rr_tol = float(args.bc_rr_tolerance)
oc_rr_tol = float(args.oc_rr_tolerance)

bc_probs_tol = float(args.bc_probs_tolerance)
oc_probs_tol = float(args.oc_probs_tolerance)

print("=============================================")
print("BC Risk Tolerance "+str(bc_rr_tol))
print("OC Risk Tolerance "+str(oc_rr_tol))

print("BC Probs Tolerance "+str(bc_probs_tol))
print("OC Probs Tolerance "+str(oc_probs_tol))


irates = "BOADICEA-Model-V6/Data/incidences_"+args.cancer_rates+".nml"
print('Cancer Incidence Rates: '+args.cancer_rates)
print("=============================================")

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
elif os.path.isdir(bwa):
    bwalist = [os.path.join(bwa, f) for f in os.listdir(bwa) if os.path.isfile(os.path.join(bwa, f))]
else:
    bwalist = glob_re(bwa)

# loop over canrisk files and compare results from webservices with those from the batch script
exact_matches = 0
diffs = []
for bwa in bwalist:
    try:
        cwd = tempfile.mkdtemp(prefix="canrisk_batch_")

        # run webservice
        args.tab = os.path.join(cwd, 'webservice.tab')
        runws(args, {"user_id": "end_user_id"}, bwa, ['boadicea', 'ovarian'], token, url)

        # create pedigree csv file for batch script
        rfsnames, rfs, ashkn = get_rfs(bwa)
        csvfile = os.path.join(cwd, "ped.csv")
        convert2csv(bwa, csvfile, rfsnames, rfs)

        # run batch script
        BC_BATCH_RISKS = os.path.join(cwd, "batch_boadicea_risks.out")
        BC_BATCH_PROBS = os.path.join(cwd, "batch_boadicea_probs.out")
        outs, errs = run_batch(FORTRAN, cwd, csvfile, BC_BATCH_RISKS, irates, ashkn=ashkn)
        outs, errs = run_batch(FORTRAN, cwd, csvfile, BC_BATCH_PROBS, irates, ashkn=ashkn, muts=True)
        OC_BATCH_RISKS = os.path.join(cwd, "batch_ovarian_risks.out")
        OC_BATCH_PROBS = os.path.join(cwd, "batch_ovarian_probs.out")
        outs, errs = run_batch(FORTRAN, cwd, csvfile, OC_BATCH_RISKS, irates, ashkn=ashkn, model='OC')
        outs, errs = run_batch(FORTRAN, cwd, csvfile, OC_BATCH_PROBS, irates, ashkn=ashkn, model='OC', muts=True)

        # get results
        c_ages = get_censoring_ages(bwa)
        bc_ws, oc_ws, bc_mp_ws, oc_mp_ws = get_ws_results(args, c_ages)

        bc_batch = get_batch_results(BC_BATCH_RISKS, c_ages)
        oc_batch = get_batch_results(OC_BATCH_RISKS, c_ages)
        bc_mp_batch = get_mp(BC_BATCH_PROBS)
        oc_mp_batch = get_mp(OC_BATCH_PROBS)

        print(bwa)
        if bc_mp_batch is not None and bc_mp_ws is not None:
            exact_matches = compare_mp("BC", bc_mp_batch, bc_mp_ws, exact_matches, bc_probs_tol)
        if oc_mp_batch is not None and oc_mp_ws is not None:
            exact_matches = compare_mp("OC", oc_mp_batch, oc_mp_ws, exact_matches, oc_probs_tol)

        # compare webservice.tab with batch_result.out
        for age in c_ages:
            if len(bc_ws) > 0 or bc_batch is not None:
                if bc_ws[age] and bc_batch[age] and math.isclose(float(bc_ws[age]), float(bc_batch[age]),
                                                                 abs_tol=bc_rr_tol):
                    print(f"BC EXACT MATCH ::: {age}    webservice: {bc_ws[age]} batch: {bc_batch[age]}", end='\t\t')
                else:
                    print(f"BC DIFFERENCE ["+str(float(bc_ws[age])-float(bc_batch[age])) +
                          f"]*** {age}    webservice: {bc_ws[age]} batch: {bc_batch[age]}", end='\t\t')
                    exact_matches += 1
                    diffs.append(bwa)

            if len(oc_ws) > 0 or oc_batch is not None:
                if oc_ws[age] and oc_batch[age] and math.isclose(float(oc_ws[age]), float(oc_batch[age]),
                                                                 abs_tol=oc_rr_tol):
                    print(f"OC EXACT MATCH :::    webservice: {oc_ws[age]} batch: {oc_batch[age]}")
                else:
                    print(f"OC DIFFERENCE ["+str(float(oc_ws[age])-float(oc_batch[age])) +
                          f"]***    webservice: {oc_ws[age]} batch: {oc_batch[age]}")
                    exact_matches += 1
                    diffs.append(bwa)
    finally:
        shutil.rmtree(cwd)
        # print(cwd)

if exact_matches != 0:
    print(f"====== DIFFERENCES FOUND {exact_matches}")
    diffs = list(set(diffs))
    for d in diffs:
        print(d)
else:
    print("====== NO DIFFERENCES FOUND")

exit(exact_matches)
