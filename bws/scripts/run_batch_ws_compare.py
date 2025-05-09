#!/usr/bin/env python3
#
# Test script to compare results from BC webservice and batch script
# bws/scripts/run_batch_ws_compare.py -u USERNAME -p bws/tests/data/canrisk_prs.txt
# bws/scripts/run_batch_ws_compare.py -u USERNAME -p bws/tests/data/batch/ --url http://0.0.0.0:8000/ \
#                                     --fortran /home/xxxx/Model-Batch-Processing --cancer_rates Spain
#
# © 2023 University of Cambridge
# SPDX-FileCopyrightText: 2023 University of Cambridge
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
from bws.risk_factors.ethnicity import UKBioBankEthnicty


def get_ws_results(args, calc_ages):
    ''' Parse web-service tab file and return breast and ovarian cancer risks and mutation carrier probabilities. '''
    bc_ws, oc_ws, pc_ws, bc_mp_ws, oc_mp_ws, pc_mp_ws = {}, {}, {}, None, None, None
    sages = '|'.join([str(c) for c in calc_ages])
    with open(args.tab, 'r') as f:
        model = 'BC'
        for line in f:
            if 'BOADICEA MODEL' in line:
                model = 'BC'
            elif 'OVARIAN MODEL' in line:
                model = 'OC'
            elif 'PROSTATE MODEL' in line:
                model = 'PC'
            crisks = re.match("^(.*\t.+\t("+sages+")\t(\d*\.?\d*[e\-\d*]*)).*", line)
            if crisks:
                if model == 'BC':
                    bc_ws[int(crisks.group(2))] = crisks.group(3)
                elif model == 'OC':
                    oc_ws[int(crisks.group(2))] = crisks.group(3)
                elif model == 'PC':
                    pc_ws[int(crisks.group(2))] = crisks.group(3)
            elif 'no mutation' in line:
                gkeys = line.strip().split('\t')[1:]    # ignore FamID in first column
                vals = next(f).strip().split('\t')[1:]
                if model == 'BC':
                    bc_mp_ws = {gkeys[idx].strip(): val for idx, val in enumerate(vals)}
                elif model == 'OC':
                    oc_mp_ws = {gkeys[idx].strip(): val for idx, val in enumerate(vals)}
                elif model == 'PC':
                    pc_mp_ws = {gkeys[idx].strip(): val for idx, val in enumerate(vals)}
        f.close()
    return bc_ws, oc_ws, pc_ws, bc_mp_ws, oc_mp_ws, pc_mp_ws


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
parser.add_argument('--mut_freq', default='UK', choices=['UK', 'UK, non-European', 'Ashkenazi', 'Iceland'],
                    help='Mutation Frequencies (default: %(default)s)')
parser.add_argument('--cancer_rates', default='UK',
                    choices=['UK', 'Australia', 'Canada', 'USA', 'Denmark', 'Estonia', 'Finland', 'France',
                             'Iceland', 'Netherlands', 'New-Zealand', 'Norway', 'Slovenia', 'Spain', 'Sweden'],
                    help='Cancer incidence rates (default: %(default)s)')
parser.add_argument('--token', help='authentication token')
parser.add_argument('--bc_rr_tolerance', default=1e-09, help='BC tolerance comparing web-service & batch risks')
parser.add_argument('--oc_rr_tolerance', default=1e-09, help='OC tolerance comparing web-service & batch risks')
parser.add_argument('--pc_rr_tolerance', default=1e-09, help='PC tolerance comparing web-service & batch risks')

parser.add_argument('--bc_probs_tolerance', default=1e-09, help='BC tolerance comparing web-service & batch probs')
parser.add_argument('--oc_probs_tolerance', default=1e-09, help='OC tolerance comparing web-service & batch probs')
parser.add_argument('--pc_probs_tolerance', default=1e-09, help='PC tolerance comparing web-service & batch probs')

args = parser.parse_args()

bc_rr_tol = float(args.bc_rr_tolerance)
oc_rr_tol = float(args.oc_rr_tolerance)
pc_rr_tol = float(args.oc_rr_tolerance)

bc_probs_tol = float(args.bc_probs_tolerance)
oc_probs_tol = float(args.oc_probs_tolerance)
pc_probs_tol = float(args.oc_probs_tolerance)

print("=============================================")
print("BC Risk Tolerance "+str(bc_rr_tol))
print("OC Risk Tolerance "+str(oc_rr_tol))
print("PC Risk Tolerance "+str(pc_rr_tol))

print("BC Probs Tolerance "+str(bc_probs_tol))
print("OC Probs Tolerance "+str(oc_probs_tol))
print("PC Probs Tolerance "+str(pc_probs_tol))

irates = args.cancer_rates.replace('New-Zealand', 'New_Zealand')+".nml"
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
mut_freq = args.mut_freq
for bwa in bwalist:
    try:
        cwd = tempfile.mkdtemp(prefix="canrisk_batch_")
        rfsnames, rfs, ashkn, biobank_ethnicity = get_rfs(bwa)

        # run webservice
        if UKBioBankEthnicty.GROUPS[biobank_ethnicity.ethnicity] ==  "UK-pop":
            oc_mut_freq = mut_freq.upper()
            bc_mut_freq = mut_freq.upper()
            pc_mut_freq = mut_freq.upper()
        else:
            oc_mut_freq = mut_freq.upper()
            bc_mut_freq = 'UK, non-European'
            pc_mut_freq = mut_freq.upper()

        args.mut_freq = bc_mut_freq
        args.tab = os.path.join(cwd, 'webservice.tab')
        runws(args, {"user_id": "end_user_id"}, bwa, ['boadicea', 'ovarian', 'prostate'], token, url)

        # create pedigree csv file for batch script
        csvfile = os.path.join(cwd, "ped.csv")
        convert2csv(bwa, csvfile, rfsnames, rfs)

        # run batch script
        BC_BATCH_RISKS = os.path.join(cwd, "batch_boadicea_risks.out")
        BC_BATCH_PROBS = os.path.join(cwd, "batch_boadicea_probs.out")
        outs, errs = run_batch(FORTRAN, cwd, csvfile, BC_BATCH_PROBS, irates, ashkn=ashkn, mut_freq=bc_mut_freq, biobank_ethnicity=biobank_ethnicity, muts=True)
        outs, errs = run_batch(FORTRAN, cwd, csvfile, BC_BATCH_RISKS, irates, ashkn=ashkn, mut_freq=bc_mut_freq, biobank_ethnicity=biobank_ethnicity)
        OC_BATCH_RISKS = os.path.join(cwd, "batch_ovarian_risks.out")
        OC_BATCH_PROBS = os.path.join(cwd, "batch_ovarian_probs.out")
        outs, errs = run_batch(FORTRAN, cwd, csvfile, OC_BATCH_RISKS, irates, ashkn=ashkn, mut_freq=oc_mut_freq, model='OC', biobank_ethnicity=biobank_ethnicity)
        outs, errs = run_batch(FORTRAN, cwd, csvfile, OC_BATCH_PROBS, irates, ashkn=ashkn, mut_freq=oc_mut_freq, model='OC', biobank_ethnicity=biobank_ethnicity, muts=True)
        PC_BATCH_RISKS = os.path.join(cwd, "batch_prostate_risks.out")
        outs, errs = run_batch(FORTRAN, cwd, csvfile, PC_BATCH_RISKS, irates, ashkn=ashkn, mut_freq=pc_mut_freq, model='PC', biobank_ethnicity=biobank_ethnicity)

        # get results
        c_ages = get_censoring_ages(bwa)
        bc_ws, oc_ws, pc_ws, bc_mp_ws, oc_mp_ws, pc_mp_ws = get_ws_results(args, c_ages)

        bc_batch = get_batch_results(BC_BATCH_RISKS, c_ages)
        oc_batch = get_batch_results(OC_BATCH_RISKS, c_ages)
        pc_batch = get_batch_results(PC_BATCH_RISKS, c_ages)
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

            if len(pc_ws) > 0 or pc_batch is not None:
                if pc_ws[age] and pc_batch[age] and math.isclose(float(pc_ws[age]), float(pc_batch[age]),
                                                                 abs_tol=pc_rr_tol):
                    print(f"PC EXACT MATCH ::: {age}    webservice: {pc_ws[age]} batch: {pc_batch[age]}")
                else:
                    print(f"PC DIFFERENCE ["+str(float(pc_ws[age])-float(pc_batch[age])) +
                          f"]*** {age}    webservice: {pc_ws[age]} batch: {pc_batch[age]}")
                    exact_matches += 1
                    diffs.append(bwa)
        if exact_matches != 0:
            break
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
