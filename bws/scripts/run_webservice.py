#!/usr/bin/env python3
#
# Requirements: requests module
#    https://pypi.org/project/requests/
#
# Run risk predictions via the web-service.
#
# optional arguments:
#  -h, --help            show this help message and exit
#  --mut_freq {UK,Ashkenazi,Iceland}
#                        Mutation Frequencies (default: UK)
#  --cancer_rates {'UK','Australia','Canada','USA','Denmark','Estonia','Finland','France',
#                  'Iceland','Netherlands','New-Zealand','Norway','Slovenia','Spain','Sweden'}
#                        Cancer incidence rates (default: UK)
#  --url URL             Web-services URL
#  -u USER, --user USER  Username
#  -p PEDIGREE, --pedigree PEDIGREE
#                        Pedigree file
#
# Genetic test sensitivity:
#  --brca1_mut_sensitivity BRCA1_MUT_SENSITIVITY
#                        brca1 mutation sensitivity
#  --brca2_mut_sensitivity BRCA2_MUT_SENSITIVITY
#                        brca2 mutation sensitivity
#  --palb2_mut_sensitivity PALB2_MUT_SENSITIVITY
#                        palb2 mutation sensitivity
#  --chek2_mut_sensitivity CHEK2_MUT_SENSITIVITY
#                        chek2 mutation sensitivity
#  --atm_mut_sensitivity ATM_MUT_SENSITIVITY
#                        atm mutation sensitivity
#
# e.g.
# run_webservice.py -u username -p ~/bwa4_beta_pedigree_data.txt
#
# run_webservice.py -u username -p bws/bws/tests/data/canrisk_multi4x.txt -c both -t example.tab
#
# run_webservice.py -u username -p boadicea/tests_selenium/canrisk_format_data/canrisk_data1.txt \
#       --vcf sample_data/sample_BCAC_313.vcf -s SampleA --bc_prs_reference_file BCAC_313_PRS.prs
#

import getpass
import json
import argparse
import csv
import os
import sys
from pathlib import Path
import pdf_report

try:
    import grequests
except ImportError as e:
    import requests

from os import listdir
from os.path import join, isfile


bc_genes = ['brca1', 'brca2', 'palb2', 'chek2', 'atm', 'bard1']
oc_genes = ['brca1', 'brca2', 'rad51c', 'rad51d', 'brip1']
pc_genes = ['hoxb13']

def post_requests(url, **kwargs):
    ''' Post requests via <code>grequests</code> if installed otherwise via <code>requests</code>. '''
    if 'grequests' in sys.modules:
        r = grequests.post(url, **kwargs)
        return grequests.map([r])[0]
    return requests.post(url, **kwargs)


def mutation_probability_output(res, writer):
    ''' Write out mutation carrier probabilities '''
    if 'mutation_probabilties' in res:
        writer.writerow(['Mutation Carrier Probabilties'])
        muts = res['mutation_probabilties']
        col = ['FamID']
        val = [res['family_id']]
        for m in muts:
            k = list(m.keys())[0]
            col.append(k+"    "+("  " if k == "ATM" else ""))
            val.append(m[k]['decimal'])
        writer.writerow(col)
        writer.writerow(val)
        writer.writerow([])


def summary_output_tab(tabf, cmodel, rjson, bwa):
    ''' Tab delimited output file '''
    if 'version' not in rjson:
        return
    if cmodel == "boadicea":
        header = ["FamID", "IndivID", "Age",
                  "+5 BC Risk", "+10 BC Risk", "80 BC Risk", "BC Lifetime"]
        ctype = "breast cancer risk"
    elif cmodel == "ovarian":
        header = ["FamID", "IndivID", "Age",
                  "+5 OC Risk", "+10 OC Risk", "80 OC Risk", "OC Lifetime"]
        ctype = "ovarian cancer risk"
    elif cmodel == "prostate":
        header = ["FamID", "IndivID", "Age",
                  "+5 PC Risk", "+10 PC Risk", "80 PC Risk", "PC Lifetime"]
        ctype = "prostate cancer risk"

    with open(tabf, 'a') as csvfile:
        writer = csv.writer(csvfile, delimiter='\t')
        writer.writerow([rjson["version"].upper()])
        writer.writerow(["===================================="])
        writer.writerow(["Pedigree: "+bwa])
        writer.writerow([rjson["timestamp"]])
        writer.writerow(["Cancer incidence rates: "+rjson["cancer_incidence_rates"]])
        writer.writerow(["Pathogenic Variant Frequency: "+list(rjson["mutation_frequency"].keys())[0]])

        results = rjson["pedigree_result"]
        for res in results:
            famid = res["family_id"]
            indivID = res["proband_id"]

            age = "-"
            c5, c10, c80, clt = "-", "-", "-", "-"

            if "cancer_risks" in res:
                writer.writerow([])
                writer.writerow(header)
                cancer_risks = res["cancer_risks"]
                for cr in cancer_risks:
                    if age == "-" or int(cr["age"]) < age:
                        age = int(cr["age"]) - 1

                for cr in cancer_risks:
                    crage = int(cr["age"])
                    if crage == (age + 5):
                        c5 = cr[ctype]["decimal"]
                    elif crage == (age + 10):
                        c10 = cr[ctype]["decimal"]
                    if crage == 80:
                        c80 = cr[ctype]["decimal"]

            # report lifetime cancer risks
            if "lifetime_cancer_risk" in res and res["lifetime_cancer_risk"] is not None:
                cr = res['lifetime_cancer_risk'][0]
                clt = cr["breast cancer risk"]["decimal"]

            writer.writerow([famid, indivID, age, c5, c10, c80, clt])
            mutation_probability_output(res, writer)
        writer.writerow([])
    csvfile.close()


def output_tab(tabf, cmodel, rjson, bwa):
    ''' Tab delimited output file '''
    if 'version' not in rjson:
        return
    if cmodel == "boadicea":
        header = ["FamID", "IndivID", "Age", "BCRisk          ", "BCRisk%"]
        cname = "breast"
    elif cmodel == "ovarian":
        header = ["FamID", "IndivID", "Age", "OCRisk          ", "OCRisk%"]
        cname = "ovarian"
    elif cmodel == "prostate":
        header = ["FamID", "IndivID", "Age", "PCRisk          ", "PCRisk%"]
        cname = "prostate"
        
    with open(tabf, 'a') as csvfile:
        writer = csv.writer(csvfile, delimiter='\t')
        writer.writerow([rjson["version"].upper()])
        writer.writerow(["===================================="])
        writer.writerow(["Pedigree", bwa])
        writer.writerow([rjson["timestamp"]])
        writer.writerow(["Cancer Incidence Rates", rjson["cancer_incidence_rates"]])
        writer.writerow(["Pathogenic Variant Frequency", list(rjson["mutation_frequency"].keys())[0]])
        writer.writerow(["NOTE: baseline cancer risks are provided in brackets"])
        results = rjson["pedigree_result"]
        for res in results:
            famid = res["family_id"]
            indivID = res["proband_id"]

            if "cancer_risks" in res:
                cancer_risks = res["cancer_risks"]
                bcancer_risks = res["baseline_cancer_risks"]

                writer.writerow(header)
                for idx, cr in enumerate(cancer_risks):
                    bcr = bcancer_risks[idx]
                    if bcr["age"] == cr["age"]:
                        cr_dec = '{} ({})'.format(cr[cname+" cancer risk"]["decimal"],
                                                  bcr[cname+" cancer risk"]["decimal"])
                        cr_per = '{} ({})'.format(cr[cname+" cancer risk"]["percent"],
                                                  bcr[cname+" cancer risk"]["percent"])
                        writer.writerow([famid, indivID, cr["age"], cr_dec, cr_per])

            # report lifetime cancer risks
            if "lifetime_cancer_risk" in res and res["lifetime_cancer_risk"] is not None:
                cr = res['lifetime_cancer_risk'][0]
                bcr = res['baseline_lifetime_cancer_risk'][0]

                cr_dec = '{} ({})'.format(cr[cname+" cancer risk"]["decimal"],
                                          bcr[cname+" cancer risk"]["decimal"])
                cr_per = '{} ({})'.format(cr[cname+" cancer risk"]["percent"],
                                          bcr[cname+" cancer risk"]["percent"])
                writer.writerow([cname.title()+" Cancer Lifetime Risk", cr_dec, cr_per])

            writer.writerow([])

            if cmodel != "prostate":
                mutation_probability_output(res, writer)

    csvfile.close()


def runws(args, data, bwa, cancers, token, url, cwd=os.getcwd(), prs=None):
    ''' Call web-services '''
    bwa = join(cwd, bwa)
    data["mut_freq"] = args.mut_freq
    data["cancer_rates"] = args.cancer_rates

    if 'grequests' in sys.modules:
        print("ASYNC")
        reqs = []
        for cmodel in cancers:
            files = {'pedigree_data': open(bwa, 'rb')}

            if 'prs' in data:
                del data['prs']
            if prs is not None:
                if cmodel == 'boadicea' and 'breast_cancer_prs' in prs and prs['breast_cancer_prs']['alpha'] != 0.0:
                    data['prs'] = json.dumps(prs['breast_cancer_prs'])
                elif cmodel == 'ovarian' and 'ovarian_cancer_prs' in prs and prs['ovarian_cancer_prs']['alpha'] != 0.0:
                    data['prs'] = json.dumps(prs['ovarian_cancer_prs'])

            reqs.append(grequests.post(url+cmodel+'/', data=data.copy(), files=files,
                                       headers={'Authorization': "Token "+token}))
        res = grequests.map(reqs)

        for idx, cmodel in enumerate(cancers):
            handle_response(args, cmodel, res[idx], bwa)
        return
    else:
        print("SYNC")
        combine = {}
        for cmodel in cancers:
            files = {'pedigree_data': open(bwa, 'rb')}

            if 'prs' in data:
                del data['prs']
            if prs is not None:
                if cmodel == 'boadicea' and 'breast_cancer_prs' in prs and prs['breast_cancer_prs']['alpha'] != 0.0:
                    data['prs'] = json.dumps(prs['breast_cancer_prs'])
                elif cmodel == 'ovarian' and 'ovarian_cancer_prs' in prs and prs['ovarian_cancer_prs']['alpha'] != 0.0:
                    data['prs'] = json.dumps(prs['ovarian_cancer_prs'])

            r = requests.post(url+cmodel+'/', data=data, files=files, headers={'Authorization': "Token "+token})
            handle_response(args, cmodel, r, bwa)
            combine[cmodel] = r
        if 'pdf' in args and args.pdf:
            pdf_report.create_pdf(url, token, combine['ovarian'], combine['boadicea'], bwa, cwd)


def handle_response(args, cmodel, r, bwa):
    ''' Handle response from web-service '''
    if r.status_code == 200:
        rjson = r.json()
        if args.tab:
            output_tab(args.tab, cmodel, rjson, bwa)
        elif args.summary:
            summary_output_tab(args.summary, cmodel, rjson, bwa)
        elif not ('pdf' in args and args.pdf):
            print(json.dumps(rjson, indent=4, sort_keys=True))
    else:
        sys.stderr.write("Web-services error status: "+str(r.status_code))
        sys.stderr.write(r.text)
        exit(1)


def get_auth_token(args, url):
    ''' Returns an authentication token. '''
    if args.token is None:
        if args.user is None:
            user = input("Username: ")
        else:
            user = args.user
        pwd = getpass.getpass()

        r = post_requests(url+"auth-token/", data={"username": user, "password": pwd})
        if r.status_code == 200:
            if args.showtoken:
                print(r.json()['token'])
                exit(0)
            return r.json()['token']
        else:
            print("Authentication error status: "+str(r.status_code))
            print(r.text)
            exit(1)
    else:
        return args.token


if __name__ == "__main__":
    cwd = Path.cwd()
    #
    # define optional command line arguments
    parser = argparse.ArgumentParser('run a risk prediction via the web-service')
    parser.add_argument('-c', '--can', default='boadicea', choices=['boadicea', 'ovarian', 'prostate', 'all'],
                        help='Cancer risk models')

    # VCF to PRS
    group1 = parser.add_argument_group('Polygenic Risk Score (PRS)')
    group1.add_argument('-v', '--vcf', help='Variant Call Format (VCF) file')
    group1.add_argument('-s', '--sample', help='sample name')
    group1.add_argument('--bc_prs_reference_file', help='breast cancer prs reference file', default=None,
                        choices=[None, 'BCAC_313_PRS.prs', 'PERSPECTIVE_295_PRS.prs', 'BRIDGES_306_PRS.prs'])
    group1.add_argument('--oc_prs_reference_file', help='ovarian cancer prs reference file', default=None,
                        choices=[None, 'OCAC_36_PRS.prs'])
    group1.add_argument('--vcfonly', help='Only run VCF to PRS', action='store_true')

    # Mutation frequencies
    parser.add_argument('--mut_freq', default='UK', choices=['UK', 'UK, non-European', 'Ashkenazi', 'Iceland'],
                        help='Mutation Frequencies (default: %(default)s)')

    genes = list(set(bc_genes + oc_genes + pc_genes))

    group2 = parser.add_argument_group('Genetic test sensitivity')
    for gene in genes:
        group2.add_argument('--'+gene+'_mut_sensitivity', type=float, help=gene+' mutation sensitivity')

    parser.add_argument('--cancer_rates', default='UK',
                        choices=['UK', 'Australia', 'Canada', 'USA', 'Denmark', 'Estonia', 'Finland', 'France',
                                 'Iceland', 'Netherlands', 'New-Zealand', 'Norway', 'Slovenia', 'Spain', 'Sweden'],
                        help='Cancer incidence rates (default: %(default)s)')

    parser.add_argument('--url', default='https://canrisk.org/', help='Web-services URL')
    parser.add_argument('-u', '--user', help='Username')
    parser.add_argument('-p', '--ped', help='CanRisk (or BOADICEA v4) pedigree file or directory of pedigree file(s)')
    parser.add_argument('-t', '--tab', help='Tab delimited output file name')
    parser.add_argument('--summary', help='Tab delimited summary output file name')
    parser.add_argument('--pdf', help='PDF file name [EXPERIMENTAL OPTION - USE WITH CAUTION]', action='store_true')
    parser.add_argument('--token', help='authentication token')
    parser.add_argument('--showtoken', help='display the authentication token', action='store_true')

    #######################################################
    args = parser.parse_args()
    if args.tab or args.summary:
        tabout = Path(args.tab if args.tab else args.summary)
        if tabout.exists():
            print("Output file already exists!")
            exit(1)

    if args.can == "all":
        cancers = ['boadicea', 'ovarian', 'prostate']
        genes = list(set(bc_genes + oc_genes + pc_genes))
    elif args.can == "ovarian":
        cancers = ['ovarian']
        genes = oc_genes
    elif args.can == "boadicea":
        cancers = ['boadicea']
        genes = bc_genes
    elif args.can == "prostate":
        cancers = ['prostate']
        genes = pc_genes

    data = {"user_id": "end_user_id"}

    for gene in genes:
        if args.__dict__[gene+"_mut_sensitivity"] is not None:
            data[gene+"_mut_sensitivity"] = args.__dict__[gene+"_mut_sensitivity"]

    url = args.url

    #######################################################
    # 1. request an authentication token

    token = get_auth_token(args, url)

    #######################################################
    # 2. optionally get PRS from VCF
    prs = None
    if args.vcf is not None:
        prs_data = {}
        if args.sample is None:
            parser.error('The vcf argument requires the sample name option [-s or --sample]')
        prs_data['sample_name'] = args.sample
        if args.bc_prs_reference_file is None and args.oc_prs_reference_file is None:
            parser.error('The vcf argument requires the breast and/or ovarian PRS reference ' +
                         'file option [--bc_prs_reference_file, --oc_prs_reference_file]')

        if args.bc_prs_reference_file is not None:
            prs_data['bc_prs_reference_file'] = args.bc_prs_reference_file
        if args.oc_prs_reference_file is not None:
            prs_data['oc_prs_reference_file'] = args.oc_prs_reference_file

        files = {'vcf_file': open(args.vcf, 'rb')}
        r = post_requests(url+'vcf2prs/', data=prs_data, files=files, headers={'Authorization': "Token "+token})
        if r.status_code == 200:
            prs = r.json()
        if args.vcfonly is not None:
            print(prs)
            exit(0)

    #######################################################
    # 3. run risk prediction for a given pedigree or all files in a directory
    bwa = (input("Pedigree (BOADICEA v4/CanRisk file) or path to directory of pedigrees: ")
           if args.ped is None else args.ped)

    # if bwa is a directory then get files in directory
    bwas = [join(bwa, f) for f in listdir(bwa) if isfile(join(bwa, f))] if os.path.isdir(bwa) else [bwa]

    if prs is not None and len(bwas) > 1:
        print("The --vcf option generates a PRS code that is expected to be used with a single pedigree file.")
        exit(1)

    http_server = None
    try:
        if 'pdf' in args and args.pdf:
            http_server = pdf_report.HttpServer()
            http_server.start_www(url)
        
        for bwa in bwas:
            print("PROCESSING: "+bwa)
            runws(args, data, bwa, cancers, token, url, cwd=cwd, prs=prs)
    finally:
        if 'pdf' in args and args.pdf:
            http_server.stop_www()
