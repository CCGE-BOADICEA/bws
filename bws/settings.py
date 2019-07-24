from collections import OrderedDict
import os

# FORTRAN settings
FORTRAN_HOME = "/home/MINTS/tjc29/boadicea_classic/"
FORTRAN_TIMEOUT = 60*4   # seconds
CWD_DIR = "/tmp"

# Environment variables for OpenBLAS (http://www.openblas.net)
FORTRAN_ENV = os.environ.copy()
FORTRAN_ENV['LD_LIBRARY_PATH'] = (FORTRAN_ENV['LD_LIBRARY_PATH']
                                  if 'LD_LIBRARY_PATH' in FORTRAN_ENV else "") + ":/usr/local/lib"
FORTRAN_ENV['OMP_STACKSIZE'] = '2G'
# FORTRAN_ENV['OPENBLAS_NUM_THREADS'] = '1'

# wkhtmltopdf executable used to generate PDF from HTML
# WKHTMLTOPDF = '/usr/bin/wkhtmltopdf'
# WKHTMLTOPDF_TIMEOUT = 10  # seconds

MAX_PEDIGREE_SIZE = 275
MIN_BASELINE_PEDIGREE_SIZE = 1
MENDEL_NULL_YEAR_OF_BIRTH = 9999

# Maximum age for risk calculation
MAX_AGE_FOR_RISK_CALCS = 79

MIN_YEAR_OF_BIRTH = 1850
BOADICEA_PEDIGREE_FORMAT_FOUR_DATA_FIELDS = 32
BOADICEA_CANRISK_FORMAT_ONE_DATA_FIELDS = 26
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
    'PROBS_EXE': 'boadicea_probs.exe',
    'RISKS_EXE': 'boadicea_risks.exe',
    'CANCERS': ['bc1', 'bc2', 'oc', 'prc', 'pac'],          # NOTE: order used by fortran pedigree file
    'GENES': ['BRCA1', 'BRCA2', 'PALB2', 'CHEK2', 'ATM'],   # NOTE: order used by fortran pedigree file
    'CALCS': ALLOWED_CALCS,
    'MUTATION_FREQUENCIES': OrderedDict([(
        'UK', {
            'BRCA1': 0.0006394,
            'BRCA2': 0.00102,
            'PALB2': 0.000575,
            'ATM': 0.001921,
            'CHEK2': 0.002614
        }),
        ('Ashkenazi', {
            'BRCA1': 0.008,
            'BRCA2': 0.006,
            'PALB2': 0.000575,
            'ATM': 0.001921,
            'CHEK2': 0.002614
        }),
        ('Iceland', {
            'BRCA1': 0.0006394,
            'BRCA2': 0.003,
            'PALB2': 0.000575,
            'ATM': 0.001921,
            'CHEK2': 0.002614
        }),
        ('Custom', {})
    ]),
    # Default genetic test sensitivities; updated BRCA1/2 as agreed with AA 10/04/17
    'GENETIC_TEST_SENSITIVITY': {
        "BRCA1": 0.9,
        "BRCA2": 0.9,
        "PALB2": 0.9,
        "ATM": 0.9,
        "CHEK2": 1.0,
    },
    # cancer incidence rate display name and corresponding file name
    'CANCER_RATES': OrderedDict([
        ('UK', 'UK'),
        # ('UK-version-1', 'UKold'),
        ('Australia', 'Australia'),
        ('Canada', 'Canada'),
        ('USA', 'USA'),
        ('Denmark', 'Denmark'),
        ('Finland', 'Finland'),
        ('Iceland', 'Iceland'),
        ('New-Zealand', 'NewZealand'),
        ('Norway', 'Norway'),
        ('Spain', 'Spain'),
        ('Sweden', 'Sweden'),
        ('Other', 'UK')
    ]),
    'PRS_REFERENCE_FILES': OrderedDict([
        ('BCAC 313', 'BCAC_313_PRS.prs'),
        ('BRIDGES 306', 'BRIDGES_306_PRS.prs'),
        ('PERSPECTIVE 295', 'PERSPECTIVE_295_PRS.prs')
    ])
}

#
# OVARIAN CANCER MODEL
OC_MODEL = {
    'NAME': 'OC',
    'HOME': os.path.join(FORTRAN_HOME, 'ovarian'),
    'PROBS_EXE': 'ovarian_probs.exe',
    'RISKS_EXE': 'ovarian_risks.exe',
    'CANCERS': ['bc1', 'bc2', 'oc', 'prc', 'pac'],              # NOTE: order used by fortran pedigree file
    'GENES': ['BRCA1', 'BRCA2', 'RAD51D', 'RAD51C', 'BRIP1'],   # NOTE: order used by fortran pedigree file
    'CALCS': ['carrier_probs', 'remaining_lifetime'],
    'MUTATION_FREQUENCIES': OrderedDict([(
        'UK', {
            'BRCA1': 0.0007947,
            'BRCA2': 0.002576,
            'RAD51D': 0.00018,
            'RAD51C': 0.00036,
            'BRIP1': 0.00044
        }),
        ('Ashkenazi', {
            'BRCA1': 0.008,
            'BRCA2': 0.006,
            'RAD51D': 0.00018,
            'RAD51C': 0.00036,
            'BRIP1': 0.00044
        }),
        ('Iceland', {
            'BRCA1': 0.0006394,
            'BRCA2': 0.003,
            'RAD51D': 0.00018,
            'RAD51C': 0.00036,
            'BRIP1': 0.00044
        }),
        ('Custom', {})
    ]),
    # Default genetic test sensitivities
    'GENETIC_TEST_SENSITIVITY': {
        "BRCA1": 0.9,
        "BRCA2": 0.9,
        "RAD51D": 0.9,
        "RAD51C": 0.9,
        "BRIP1": 0.9
    },
    # cancer incidence rate display name and corresponding file name
    'CANCER_RATES': OrderedDict([
        ('UK', 'UK'),
        # ('UK-version-1', 'UKold'),
        ('Australia', 'Australia'),
        ('Canada', 'Canada'),
        ('USA', 'USA'),
        ('Denmark', 'Denmark'),
        ('Finland', 'Finland'),
        ('Iceland', 'Iceland'),
        ('New-Zealand', 'NewZealand'),
        ('Norway', 'Norway'),
        ('Spain', 'Spain'),
        ('Sweden', 'Sweden'),
        ('Other', 'UK')
    ]),
    'PRS_REFERENCE_FILES': OrderedDict()
}

# Minimum allowable BRCA1/2 mutation is set to 0.0001. We should not allow zero, because if
# there is a mutation it will give an inconsistency, and we know that zero is unrealistic.
MIN_MUTATION_FREQ = 0.0001
# Maximum mutation frequency is currently set to Ashkenazi BRCA1 mutation value
MAX_MUTATION_FREQ = 0.008

# rest framework throttle settings
REST_FRAMEWORK = {
    'DEFAULT_THROTTLE_RATES': {
        'sustained': '5000/day',
        'burst': '300/min',
        'enduser_burst': '300/min'
    }
}
