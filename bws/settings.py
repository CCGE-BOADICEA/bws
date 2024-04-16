"""
© 2023 University of Cambridge
SPDX-FileCopyrightText: 2023 University of Cambridge
SPDX-License-Identifier: GPL-3.0-or-later
"""

from collections import OrderedDict
import os
import re
from pathlib import Path

try:
    import vcf2prs
except ImportError as e:
    pass  # module doesn't exist, deal with it.


def line_that_contain(s1, s2, fp):
    return [line for line in fp if s1 in line and s2 in line][0]


def get_alpha(ref_file):
    ''' Get PRS alpha from a reference file header. '''
    alpha = ""
    try:
        moduledir = Path(vcf2prs.__file__).parent.parent
        ref_file = os.path.join(moduledir, "PRS_files", ref_file)
        try:
            snp_file = open(ref_file, 'r')
            alpha = vcf2prs.prsinfo.PrsInfo.extract_alpha(line_that_contain('alpha', "=", snp_file))
        except (IOError, UnicodeDecodeError, StopIteration, vcf2prs.exception.Vcf2PrsError):
            raise vcf2prs.exception.Vcf2PrsError('Error: Unable to open the file "{0}".'.format(ref_file))
        finally:
            snp_file.close()
    except NameError:
        pass 
    return alpha


# FORTRAN settings
FORTRAN_HOME = "/home/tim/boadicea/"
FORTRAN_TIMEOUT = 60*4   # seconds
CWD_DIR = "/tmp"

# Environment variables for OpenBLAS (http://www.openblas.net)
FORTRAN_ENV = os.environ.copy()
FORTRAN_ENV['LD_LIBRARY_PATH'] = (FORTRAN_ENV['LD_LIBRARY_PATH']
                                  if 'LD_LIBRARY_PATH' in FORTRAN_ENV else "") + ":/usr/local/lib"
FORTRAN_ENV['OMP_STACKSIZE'] = '10M'
FORTRAN_ENV['OPENBLAS_NUM_THREADS'] = '1'

# wkhtmltopdf executable used to generate PDF from HTML
# WKHTMLTOPDF = '/usr/bin/wkhtmltopdf'
# WKHTMLTOPDF_TIMEOUT = 10  # seconds
REGEX_ASHKN = re.compile("^(Ashkenazi)$")

MAX_PEDIGREE_SIZE = 275
MIN_BASELINE_PEDIGREE_SIZE = 1
MENDEL_NULL_YEAR_OF_BIRTH = -1

# Maximum age for risk calculation
MAX_AGE_FOR_RISK_CALCS = 79

MIN_YEAR_OF_BIRTH = 1850
BOADICEA_PEDIGREE_FORMAT_FOUR_DATA_FIELDS = 32
BOADICEA_CANRISK_FORMAT_ONE_DATA_FIELDS = 26
BOADICEA_CANRISK_FORMAT_TWO_DATA_FIELDS = 27
MAX_LENGTH_PEDIGREE_NUMBER_STR = 13
MIN_FAMILY_ID_STR_LENGTH = 1
MAX_FAMILY_ID_STR_LENGTH = 7
MAX_AGE = 125
MAX_NUMBER_OF_SIBS_PER_NUCLEAR_FAMILY = 20
MAX_NUMBER_OF_SIBS_PER_NUCLEAR_FAMILY_WITH_SAME_YOB = 12
MAX_NUMBER_MZ_TWIN_PAIRS = 10
UNIQUE_TWIN_IDS = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "A"]

#
# allowed model calculations
ALLOWED_CALCS = ['carrier_probs', 'remaining_lifetime', "lifetime", "ten_year"]

#
# BREAST CANCER MODEL
BC_MODEL = {
    'NAME': 'BC',
    'HOME': os.path.join(FORTRAN_HOME, 'boadicea'),
    'EXE': 'boadicea.exe',
    'CANCERS': ['bc1', 'bc2', 'oc', 'prc', 'pac'],          # NOTE: order used by fortran pedigree file
    'GENES': ['BRCA1', 'BRCA2', 'PALB2', 'CHEK2',           # NOTE: order used by fortran pedigree file
              'ATM', 'BARD1', 'RAD51C', 'RAD51D'],
    'CALCS': ALLOWED_CALCS,
    'MUTATION_FREQUENCIES': OrderedDict([(
        'UK', {
            'BRCA1': 0.0006394,
            'BRCA2': 0.00102,
            'PALB2': 0.00064,
            'ATM': 0.0018,
            'CHEK2': 0.00373,
            'BARD1': 0.00043,
            'RAD51C': 0.00035,
            'RAD51D': 0.00035
        }),
        ('UK, non-European', {
            'BRCA1': 0.0006394,
            'BRCA2': 0.00102,
            'PALB2': 0.00064,
            'ATM': 0.0018,
            'CHEK2': 0.00109,
            'BARD1': 0.00043,
            'RAD51C': 0.00035,
            'RAD51D': 0.00035
        }),
        ('Ashkenazi', {
            'BRCA1': 0.008,
            'BRCA2': 0.006,
            'PALB2': 0.00064,
            'ATM': 0.0018,
            'CHEK2': 0.00373,
            'BARD1': 0.00043,
            'RAD51C': 0.00035,
            'RAD51D': 0.00035
        }),
        ('Iceland', {
            'BRCA1': 0.0006394,
            'BRCA2': 0.003,
            'PALB2': 0.00064,
            'ATM': 0.0018,
            'CHEK2': 0.00373,
            'BARD1': 0.00043,
            'RAD51C': 0.00035,
            'RAD51D': 0.00035
        })
    ]),
    # Default genetic test sensitivities; updated BRCA1/2 as agreed with AA 10/04/17
    'GENETIC_TEST_SENSITIVITY': {
        "BRCA1": 0.89,
        "BRCA2": 0.96,
        "PALB2": 0.92,
        "ATM": 0.94,
        "CHEK2": 0.98,
        "BARD1": 0.89,
        "RAD51C": 0.78,
        "RAD51D": 0.86,
    },
    # cancer incidence rate display name and corresponding file name
    'CANCER_RATES': OrderedDict([
        ('UK', 'UK'),
        # ('UK-version-1', 'UKold'),
        ('Australia', 'Australia'),
        ('Canada', 'Canada'),
        ('USA', 'USA'),
        ('Denmark', 'Denmark'),
        ('Estonia', 'Estonia'),
        ('Finland', 'Finland'),
        ('France', 'France'),
        ('Iceland', 'Iceland'),
        ('Netherlands', 'Netherlands'),
        ('New-Zealand', 'New_Zealand'),
        ('Norway', 'Norway'),
        ('Slovenia', 'Slovenia'),
        ('Spain', 'Spain'),
        ('Sweden', 'Sweden'),
        ('Other', 'UK')
    ]),
    'PRS_REFERENCE_FILES': OrderedDict([
        # ('ANTEBC 2803', {'alpha': xxxx}),
        ('BCAC 77', 'BCAC_77_PRS.prs'),
        ('BCAC 313 - European', 'BCAC_313_PRS.prs'),
        ('BCAC 3820', 'BCAC_3820_PRS.prs'),
        ('BRIDGES 306', 'BRIDGES_306_PRS.prs'),
        # ('DBDS 299', 'DBDS_299_PRS.prs'),
        ('EGLH-CEN 301', 'EGLH-CEN_301_PRS.prs'),
        ('EGLH-CEN 303', 'EGLH-CEN_303_PRS.prs'),
        # ('EMERGE 309', 'EMERGE_309_PRS.prs'),
        ('PERSPECTIVE 295', 'PERSPECTIVE_295_PRS.prs'),
        ('PRISMA 268', 'PRISM_268_PRS.prs'),
        ('WISDOM 75', 'WISDOM_75_PRS.prs'),
        ('WISDOM 128', 'WISDOM_128_PRS.prs')
    ])
}
BC_MODEL["INCIDENCE"] = os.path.join(BC_MODEL["HOME"], 'Data') + "/incidences_"
BC_MODEL['PRS_ALPHA'] = {
    k: get_alpha(v) if isinstance(v, str) else v['alpha']
    for k, v in BC_MODEL['PRS_REFERENCE_FILES'].items()
}

#
# OVARIAN CANCER MODEL
OC_MODEL = {
    'NAME': 'OC',
    'HOME': os.path.join(FORTRAN_HOME, 'ovarian'),
    'EXE': 'ovarian.exe',
    'CANCERS': ['bc1', 'bc2', 'oc', 'prc', 'pac'],      # NOTE: order used by fortran pedigree file
    'GENES': ['BRCA1', 'BRCA2', 'RAD51D', 'RAD51C',     # NOTE: order used by fortran pedigree file
              'BRIP1', 'PALB2'],
    'CALCS': ['carrier_probs', 'remaining_lifetime'],
    'MUTATION_FREQUENCIES': OrderedDict([(
        'UK', {
            'BRCA1': 0.0007947,
            'BRCA2': 0.002576,
            'RAD51D': 0.00035,
            'RAD51C': 0.00035,
            'BRIP1': 0.00071,
            'PALB2': 0.00064
        }),
        ('UK, non-European', {
            'BRCA1': 0.0007947,
            'BRCA2': 0.002576,
            'RAD51D': 0.00035,
            'RAD51C': 0.00035,
            'BRIP1': 0.00071,
            'PALB2': 0.00064
        }),
        ('Ashkenazi', {
            'BRCA1': 0.008,
            'BRCA2': 0.006,
            'RAD51D': 0.00035,
            'RAD51C': 0.00035,
            'BRIP1': 0.00071,
            'PALB2': 0.00064
        }),
        ('Iceland', {
            'BRCA1': 0.0007947,
            'BRCA2': 0.003,
            'RAD51D': 0.00035,
            'RAD51C': 0.00035,
            'BRIP1': 0.00071,
            'PALB2': 0.00064
        })
    ]),
    # Default genetic test sensitivities
    'GENETIC_TEST_SENSITIVITY': {
        "BRCA1": 0.89,
        "BRCA2": 0.96,
        "RAD51D": 0.86,
        "RAD51C": 0.78,
        "BRIP1": 0.95,
        "PALB2": 0.92
    },
    # cancer incidence rate display name and corresponding file name
    'CANCER_RATES': OrderedDict([
        ('UK', 'UK'),
        # ('UK-version-1', 'UKold'),
        ('Australia', 'Australia'),
        ('Canada', 'Canada'),
        ('USA', 'USA'),
        ('Denmark', 'Denmark'),
        ('Estonia', 'Estonia'),
        ('Finland', 'Finland'),
        ('France', 'France'),
        ('Iceland', 'Iceland'),
        ('Netherlands', 'Netherlands'),
        ('New-Zealand', 'New_Zealand'),
        ('Norway', 'Norway'),
        ('Slovenia', 'Slovenia'),
        ('Spain', 'Spain'),
        ('Sweden', 'Sweden'),
        ('Other', 'UK')
    ]),
    'PRS_REFERENCE_FILES': OrderedDict([
        ('OC-EGLH-CEN 34', 'OC_EGLH-CEN_34_PRS.prs'),
        ('OCAC 36', 'OCAC_36_PRS.prs')
    ])
}
OC_MODEL["INCIDENCE"] = os.path.join(OC_MODEL["HOME"], 'Data') + "/incidences_"
OC_MODEL['PRS_ALPHA'] = {key: get_alpha(value) for key, value in OC_MODEL['PRS_REFERENCE_FILES'].items()}


#
# PROSTATE CANCER MODEL
PC_MODEL = {
    'NAME': 'PC',
    'HOME': os.path.join(FORTRAN_HOME, 'prostate'),
    'EXE': 'prostate.exe',
    'CANCERS': ['prc', 'bc1', 'oc'],        # NOTE: order used by fortran pedigree file
    'GENES': ['BRCA2', 'HOXB13', 'BRCA1'],  # NOTE: order used by fortran pedigree file
    'CALCS': ['carrier_probs', 'remaining_lifetime'],
    'MUTATION_FREQUENCIES': OrderedDict([(
        'UK', {
            'BRCA1': 0.0006394,
            'BRCA2': 0.00102,
            'HOXB13': 0.00212
        }),
        ('Ashkenazi', {
            'BRCA1': 0.008,
            'BRCA2': 0.006,
            'HOXB13': 0.00212
        }),
        ('Iceland', {
            'BRCA1': 0.0006394,
            'BRCA2': 0.003,
            'HOXB13': 0.00212
        })
    ]),
    # Default genetic test sensitivities
    'GENETIC_TEST_SENSITIVITY': {
        "BRCA1": 0.89,
        "BRCA2": 0.96,
        "HOXB13": 0.86
    },
    # cancer incidence rate display name and corresponding file name
    'CANCER_RATES': OrderedDict([
        ('UK', 'UK'),
        # ('UK-version-1', 'UKold'),
        ('Australia', 'Australia'),
        ('Canada', 'Canada'),
        ('USA', 'USA'),
        ('Denmark', 'Denmark'),
        ('Estonia', 'Estonia'),
        ('Finland', 'Finland'),
        ('France', 'France'),
        ('Iceland', 'Iceland'),
        ('Netherlands', 'Netherlands'),
        ('New-Zealand', 'New_Zealand'),
        ('Norway', 'Norway'),
        ('Slovenia', 'Slovenia'),
        ('Spain', 'Spain'),
        ('Sweden', 'Sweden'),
        ('Other', 'UK')
    ]),
    'PRS_REFERENCE_FILES': OrderedDict([
        ('OC-EGLH-CEN 34', 'OC_EGLH-CEN_34_PRS.prs'),
        ('OCAC 36', 'OCAC_36_PRS.prs')
    ])
}
PC_MODEL["INCIDENCE"] = os.path.join(PC_MODEL["HOME"], 'Data') + "/incidences_"
PC_MODEL['PRS_ALPHA'] = {key: get_alpha(value) for key, value in PC_MODEL['PRS_REFERENCE_FILES'].items()}


# Minimum allowable BRCA1/2 mutation is set to 0.0001. We should not allow zero, because if
# there is a mutation it will give an inconsistency, and we know that zero is unrealistic.
MIN_MUTATION_FREQ = 0.0001
# Maximum mutation frequency is currently set to Ashkenazi BRCA1 mutation value
MAX_MUTATION_FREQ = 0.008

# rest framework throttle settings
REST_FRAMEWORK = {
    'DEFAULT_SCHEMA_CLASS': 'rest_framework.schemas.coreapi.AutoSchema',
    'DEFAULT_THROTTLE_RATES': {
        'sustained': '6000/day',
        'burst': '250/min',
        'enduser_burst': '150/min'
    }
}
