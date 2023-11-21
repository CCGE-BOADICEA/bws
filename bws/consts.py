"""
Constants 
© 2023 University of Cambridge
SPDX-FileCopyrightText: 2023 University of Cambridge
SPDX-License-Identifier: GPL-3.0-or-later
"""
import re

# BOADICEA header
REGEX_BWA_PEDIGREE_FILE_HEADER_ONE = \
    re.compile("^(BOADICEA\\s+import\\s+pedigree\\s+file\\s+format\\s[124](.0)*)$")
REGEX_CANRISK1_PEDIGREE_FILE_HEADER = re.compile("^(##CanRisk\\s1(.0)*)$")
REGEX_CANRISK2_PEDIGREE_FILE_HEADER = re.compile("^(##CanRisk\\s2(.\\d)*)$")
REGEX_CANRISK3_PEDIGREE_FILE_HEADER = re.compile("^(##CanRisk\\s3(.0)*)$")
    
REGEX_ALPHANUM_COMMAS = re.compile("^([\\w,]+)$")
REGEX_ALPHANUM_HYPHENS = re.compile("^([\\w\-]+)$")
REGEX_ONLY_HYPHENS = re.compile("^([\-]+)$")
REGEX_ONLY_ZEROS = re.compile("^[0]+$")
REGEX_AGE = re.compile("^\\d{1,3}$")
REGEX_YEAR_OF_BIRTH = re.compile("^0|((17|18|19|20)[0-9][0-9])$")
REGEX_ASHKENAZI_STATUS = re.compile("^[01]$")

BLANK_LINE = re.compile(r'^\s*$')
