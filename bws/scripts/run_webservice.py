#!/usr/bin/env python3
#
# Run a risk prediction via the web-service.
#
# optional arguments:
#  -h, --help            show this help message and exit
#  --mut_freq {UK,Ashkenazi,Iceland,Custom}
#                        Mutation Frequencies (default: UK)
#  --cancer_rates {UK,Australia,Canada,USA-white,Denmark,Finland,Iceland,New-Zealand,Norway,Spain,Sweden}
#                        Cancer incidence rates (default: UK)
#  --url URL             Web-services URL
#  -u USER, --user USER  Username
#  -p PEDIGREE, --pedigree PEDIGREE
#                        Pedigree file
#
# Gene mutation frequencies (when --mut_freq Custom):
#  --brca1_mut_frequency BRCA1_MUT_FREQUENCY
#                        brca1 mutation frequency
#  --brca2_mut_frequency BRCA2_MUT_FREQUENCY
#                        brca2 mutation frequency
#  --palb2_mut_frequency PALB2_MUT_FREQUENCY
#                        palb2 mutation frequency
#  --chek2_mut_frequency CHEK2_MUT_FREQUENCY
#                        chek2 mutation frequency
#  --atm_mut_frequency ATM_MUT_FREQUENCY
#                        atm mutation frequency
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
# run_webservice.py --mut_freq Custom --brca1_mut_frequency 0.00064 --brca2_mut_frequency 0.00102 \
#      --palb2_mut_frequency 0.000575 --chek2_mut_frequency 0.002614 --atm_mut_frequency 0.001921
#

import getpass
import requests
import json
import argparse

#
# define optional command line arguments
parser = argparse.ArgumentParser('run a risk prediction via the web-service')
parser.add_argument('--mut_freq', default='UK', choices=['UK', 'Ashkenazi', 'Iceland', 'Custom'],
                    help='Mutation Frequencies (default: %(default)s)')

genes = ['brca1', 'brca2', 'palb2', 'chek2', 'atm']
group1 = parser.add_argument_group('Gene mutation frequencies (when --mut_freq Custom)')
for gene in genes:
    group1.add_argument('--'+gene+'_mut_frequency', type=float, help=gene+' mutation frequency')

group2 = parser.add_argument_group('Genetic test sensitivity')
for gene in genes:
    group2.add_argument('--'+gene+'_mut_sensitivity', type=float, help=gene+' mutation sensitivity')

parser.add_argument('--cancer_rates', default='UK',
                    choices=['UK', 'Australia', 'Canada', 'USA-white', 'Denmark', 'Finland',
                             'Iceland', 'New-Zealand', 'Norway', 'Spain', 'Sweden'],
                    help='Cancer incidence rates (default: %(default)s)')

parser.add_argument('--url', default='https://canrisk.org/', help='Web-services URL')
parser.add_argument('-u', '--user', help='Username')
parser.add_argument('-p', '--pedigree', help='Pedigree file')

args = parser.parse_args()
data = {"user_id": "end_user_id"}
if args.mut_freq == 'Custom':
    for gene in genes:
        if args.__dict__[gene+"_mut_frequency"] is None:
            print("--mut_freq Custom requires --"+gene+"_mut_frequency")
            exit(1)
        else:
            data[gene+"_mut_frequency"] = args.__dict__[gene+"_mut_frequency"]

for gene in genes:
    if args.__dict__[gene+"_mut_sensitivity"] is not None:
        data[gene+"_mut_sensitivity"] = args.__dict__[gene+"_mut_sensitivity"]

#
# prompt for required inputs
if args.user is None:
    user = input("Username: ")
else:
    user = args.user
pwd = getpass.getpass()
url = args.url

# 1. request an authentication token
r = requests.post(url+"auth-token/", data={"username": user, "password": pwd})
if r.status_code == 200:
    token = r.json()['token']
else:
    print("Error status: "+str(r.status_code))
    exit(1)

# 2. run BOADICEA risk prediction for a given pedigree
bwa = input("Pedigree (BOADICEA v4 file): ") if args.pedigree is None else args.pedigree

data["mut_freq"] = args.mut_freq
data["cancer_rates"] = args.cancer_rates
files = {'pedigree_data': open(bwa, 'rb')}
r = requests.post(url+"boadicea/", data=data, files=files, auth=(user, pwd))

# if successful print the results
if r.status_code == 200:
    print(json.dumps(r.json(), indent=4, sort_keys=True))
else:
    print("Error status: "+str(r.status_code))
    print(r.json())
    exit(1)
