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


def run_batch(cwd, csvfile, ofile, irates, ashkn=False, model='BC', muts=False):
    ''' Run batch processing script. '''
    setting = FORTRAN+"settings_"+model+("_AJ" if ashkn else "")+".ini"
    cmd = [FORTRAN+"batch_run.sh",
           "-r", ofile,
           "-i", irates,
           "-s", setting,
           "-l", os.path.join(cwd, model+"runlog.log")]
    if model == 'OC':
        cmd.append('-o')
    if muts:
        cmd.append('-p')
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


def get_rfs(bwa):
    '''  get risk factor names and values plus PRS '''
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
                    elif ":" in line[1]:
                        rfsnames.append(['OC_Duration', 'OC_Duration'])
                        parts = line[1].split(":")
                        rfs['OC_Duration'] = parts[1]
                        line[1] = parts[0]
                elif line[0] == 'birads':
                    name = 'BIRADS'
                elif line[0] == 'endo':
                    name = 'Endometriosis'
                else:
                    name = line[0].capitalize()

                rfsnames.append([line[0], name])
                rfs[line[0]] = line[1]
    f.close()

    with open(bwa, 'r') as f:
        pedigree_data = f.read()
    ashkn = PedigreeFile(pedigree_data).pedigrees[0].is_ashkn()
    return rfsnames, rfs, ashkn


def get_ws_results(args):
    ''' Parse web-service tab file and return breast and ovarian cancer risks and mutation carrier probabilities. '''
    bc_80_ws, oc_80_ws, bc_mp_ws, oc_mp_ws = None, None, None, None
    with open(args.tab, 'r') as f:
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
            elif line.startswith('no mutation'):
                gkeys = line.strip().split('\t')
                vals = next(f).strip().split('\t')
                if model == 'BC':
                    bc_mp_ws = {gkeys[idx]: val for idx, val in enumerate(vals)}
                else:
                    oc_mp_ws = {gkeys[idx]: val for idx, val in enumerate(vals)}
        f.close()
    return bc_80_ws, oc_80_ws, bc_mp_ws, oc_mp_ws


def get_80(fname):
    ''' Get risk at 80 from batch '''
    if not os.path.isfile(fname):
        print("BATCH FILE NOT FOUND :: "+fname)

    f = open(fname, "r")
    for line in f:
        crisks = line.strip().split(',')
        if len(crisks) == 4 and not line.startswith('file'):
            c_80_batch = crisks[3]
    f.close()
    return c_80_batch


def compare_mp(model, mp_batch, mp_ws, exact_matches):
    ''' '''
    exact = True
    msg = ""
    for k, v in mp_batch.items():
        if math.isclose(float(v), float(mp_ws[k])):
            msg += k+":"+mp_ws[k]+"="+v+" "
        else:
            msg += k+":"+mp_ws[k]+"?"+v+" "
            exact_matches += 1
            exact = False
    if exact:
        print(model+" EXACT MATCH ::: "+msg)
    else:
        print(model+" NOT A MATCH ::: "+msg)
    return exact_matches


def get_mp(fname):
    ''' Get mutation propabilities from batch '''
    if not os.path.isfile(fname):
        print("BATCH FILE NOT FOUND :: "+fname)
    if os.stat(fname).st_size == 0:
        return None

    f = open(fname, "r")
    for line in f:
        parts = line.strip().replace('NO_PATHOGENIC_VARIANTS', 'no mutation').split(',')[2:]
        if line.startswith('file'):
            keys = parts
        elif not line.startswith('file'):
            vals = parts
    f.close()
    return {keys[idx]: val for idx, val in enumerate(vals)}


# loop over canrisk files and compare results from webservices with those from the batch script
exact_matches = 0
for bwa in bwalist:
    try:
        cwd = tempfile.mkdtemp(prefix="canrisk_batch_")

        # run webservice
        args.tab = os.path.join(cwd, 'webservice.tab')
        runws(args, {"user_id": "end_user_id"}, bwa, ['boadicea', 'ovarian'], token, url)
        bc_80_ws, oc_80_ws, bc_mp_ws, oc_mp_ws = get_ws_results(args)

        # create pedigree csv file for batch script
        rfsnames, rfs, ashkn = get_rfs(bwa)
        csvfile = os.path.join(cwd, "ped.csv")
        convert2csv(bwa, csvfile, rfsnames, rfs)

        # run batch script
        BC_BATCH_RISKS = os.path.join(cwd, "batch_boadicea_risks.out")
        BC_BATCH_PROBS = os.path.join(cwd, "batch_boadicea_probs.out")
        outs, errs = run_batch(cwd, csvfile, BC_BATCH_RISKS, irates, ashkn=ashkn)
        outs, errs = run_batch(cwd, csvfile, BC_BATCH_PROBS, irates, ashkn=ashkn, muts=True)
        OC_BATCH_RISKS = os.path.join(cwd, "batch_ovarian_risks.out")
        OC_BATCH_PROBS = os.path.join(cwd, "batch_ovarian_probs.out")
        outs, errs = run_batch(cwd, csvfile, OC_BATCH_RISKS, irates, ashkn=ashkn, model='OC')
        outs, errs = run_batch(cwd, csvfile, OC_BATCH_PROBS, irates, ashkn=ashkn, model='OC', muts=True)

        bc_80_batch = get_80(BC_BATCH_RISKS)
        oc_80_batch = get_80(OC_BATCH_RISKS)
        bc_mp_batch = get_mp(BC_BATCH_PROBS)
        oc_mp_batch = get_mp(OC_BATCH_PROBS)

        if bc_mp_batch is not None or bc_mp_ws is not None:
            exact_matches = compare_mp("BC", bc_mp_batch, bc_mp_ws, exact_matches)
        if oc_mp_batch is not None or oc_mp_ws is not None:
            exact_matches = compare_mp("OC", oc_mp_batch, oc_mp_ws, exact_matches)

        # compare webservice.tab with batch_result.out
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

if exact_matches != 0:
    print("====== DIFFERENCES FOUND")
else:
    print("====== NO DIFFERENCES FOUND")

exit(exact_matches)
