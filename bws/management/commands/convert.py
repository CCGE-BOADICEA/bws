"""
Command line utility.

© 2022 Cambridge University
SPDX-FileCopyrightText: 2022 Cambridge University
SPDX-License-Identifier: GPL-3.0-or-later
"""
import os
from django.core.management.base import BaseCommand, CommandError
from bws.exceptions import PedigreeFileError
from bws.pedigree import BLANK_LINE


def convert_boadicea_pedigree_v4(pedigree_file):
    """
    Convert old BOADICEA import pedigree v1 and v2 files to v4
    @param pedigree_file: path to pedigree input file
    """
    fin = open(pedigree_file, 'r')
    fout = open(pedigree_file + '.bwa_v4.txt', 'w')
    version = ''
    for line in fin:
        if line.startswith("BOADICEA import pedigree"):
            if '2' in line:
                version = 2
            elif '1' in line:
                version = 1
            else:
                raise PedigreeFileError("Unexpected BOADICEA import pedigree version")
            fout.write(line.replace('2', '4').replace('1', '4').strip()+'\n')
        elif line.startswith("FamID"):
            fout.write("FamID\tName\tTarget\tIndivID\tFathID\tMothID\tSex\tMZtwin\tDead\tAge\tYob"
                       "\t1stBrCa\t2ndBrCa\tOvCa\tProCa\tPanCa\tAshkn\tBRCA1t\tBRCA1r\tBRCA2t\tBRCA2r"
                       "\tPALB2t\tPALB2r\tATMt\tATMr\tCHEK2t\tCHEK2r\tER\tPR\tHER2\tCK14\tCK56\n")
        elif BLANK_LINE.match(line):
            continue
        else:
            record = line.split()
            gtest = record[16]
            mutn = record[17]

            # In BOADICEA pedigree format 1, there are three possible 'Gtest' values:
            #    '0' for untested
            #    'S' for mutation search
            #    'T' for direct gene test
            #
            #    and there are five possible 'Mutn' values:
            #    '0' for no test result
            #    'N' for negative
            #    '1' for BRCA1 positive
            #    '2' for BRCA2 positive
            #    '3' for BRCA1 and BRCA2 positive
            #
            # The 15 combinations are labelled 'A' to 'O' in the table below. Table entries marked with
            # (Error) are data combinations where the input parameters 'Gtest' and 'Mutn' are inconsistent.
            # For example, in D(Error), 'Gtest' is '0' (untested) but a test result is recorded, as 'Mutn' is 'N'.
            # If a BRCA test type is specified, then the corresponding BRCA test result must be specified too
            # (and vice versa).
            #                           'Gtest' (BRCA test type)
            #                           |      0       |      S       |      T       |
            #    --------------------------------------------------------------------------------
            #                       0   |      A       |   B(Error)   |   C(Error)   |
            #     'Mutn'            N   |   D(Error)   |      E       |      F       |
            #     (BRCA test        1   |   G(Error)   |      H       |      I       |
            #     result)           2   |   J(Error)   |      K       |      L       |
            #                       3   |   M(Error)   |      N       |      O       |

            # BOADICEA format 4 BRCA genetic test results are now represented differently:
            # There are three possible 'BRCA1test'/'BRCA2test' values:
            #   '0' for untested
            #   'S' for BRCA1/BRCA2 mutation search
            #   'T' for BRCA1/BRCA2 direct gene test
            #
            #   ...and there are three possible 'BRCA1result' values:
            #   '0' for no BRCA1/BRCA2 test result
            #   'N' for BRCA1/BRCA2 negative
            #   'P' for BRCA1/BRCA2 positive
            #
            # As a result, the 15 different data combinations (labelled 'A' to 'O' in the table rows below) are now
            # represented in BOADICEA format 4 as follows
            #                    BOADICEA format 4 parameters
            #
            #                    |  BRCA1test  |  BRCA1result  |  BRCA2test  |  BRCA2result  |
            # -----------------------------------------------------------------------------------
            #                A   |      0      |       0       |      0      |       0       |
            #                B   |    Error    |       -       |      -      |       -       |
            # Data           C   |    Error    |       -       |      -      |       -       |
            # Combination    D   |    Error    |       -       |      -      |       -       |
            #                E   |      S      |       N       |      S      |       N       |
            #                F   |      T      |       N       |      T      |       N       |
            #                G   |    Error    |       -       |      -      |       -       |
            #                H   |      S      |       P       |      S      |       N       |
            #                I   |      T      |       P       |      T      |       N       |
            #                J   |    Error    |       -       |      -      |       -       |
            #                K   |      S      |       N       |      S      |       P       |
            #                L   |      T      |       N       |      T      |       P       |
            #                M   |    Error    |       -       |      -      |       -       |
            #                N   |      S      |       P       |      S      |       P       |
            #                O   |      T      |       P       |      T      |       P       |
            # -----------------------------------------------------------------------------------
            if gtest == '0' and mutn == '0':  # Test outcome A
                brca1t = '0'
                brca1r = '0'
                brca2t = '0'
                brca2r = '0'
            elif((gtest == 'S' and mutn == '0') or  # Test outcome B
                 (gtest == 'T' and mutn == '0')):  # Test outcome C
                raise PedigreeFileError("A family member has had a mutation search but no "
                                        "genetic test result has been specified")
            elif((gtest == '0' and mutn == 'N')):  # Test outcome D
                raise PedigreeFileError("A family member has been assigned a genetic test "
                                        "result but the genetic test type is unspecified")
            elif((gtest == 'S' and mutn == 'N')):  # Test outcome E
                brca1t = 'S'
                brca1r = 'N'
                brca2t = 'S'
                brca2r = 'N'
            elif((gtest == 'T' and mutn == 'N')):  # Test outcome F
                brca1t = 'T'
                brca1r = 'N'
                brca2t = 'T'
                brca2r = 'N'
            elif((gtest == '0' and mutn == '1')):  # Test outcome G
                raise PedigreeFileError("A family member has identified as BRCA1 positive "
                                        "but the genetic test type is unspecified")
            elif((gtest == 'S' and mutn == '1')):  # Test outcome H
                brca1t = 'S'
                brca1r = 'P'
                brca2t = 'S'
                brca2r = 'N'
            elif((gtest == 'T' and mutn == '1')):  # Test outcome I
                brca1t = 'T'
                brca1r = 'P'
                brca2t = 'T'
                brca2r = 'N'
            elif((gtest == '0' and mutn == '2')):  # Test outcome J
                raise PedigreeFileError("A family member has identified as BRCA2 positive "
                                        "but the genetic test type is unspecified")
            elif((gtest == 'S' and mutn == '2')):  # Test outcome K
                brca1t = 'S'
                brca1r = 'N'
                brca2t = 'S'
                brca2r = 'P'
            elif((gtest == 'T' and mutn == '2')):  # Test outcome L
                brca1t = 'T'
                brca1r = 'N'
                brca2t = 'T'
                brca2r = 'P'
            elif((gtest == '0' and mutn == '3')):  # Test outcome M
                raise PedigreeFileError("A family member has identified as BRCA1 and BRCA2 positive "
                                        "but the genetic test type is unspecified")
            elif((gtest == 'S' and mutn == '3')):  # Test outcome N
                brca1t = 'S'
                brca1r = 'P'
                brca2t = 'S'
                brca2r = 'P'
            elif((gtest == 'T' and mutn == '3')):  # Test outcome O
                brca1t = 'T'
                brca1r = 'P'
                brca2t = 'T'
                brca2r = 'P'
            else:
                raise PedigreeFileError("Program string has unexpected value")

            # Pathology results
            if version == 1:
                er, pr, her2, ck14, ck56 = '0'
            elif version == 2:
                er = record[19]
                pr = record[20]
                her2 = record[21]
                ck14 = record[22]
                ck56 = record[23]
            else:
                raise PedigreeFileError("Unexpected BOADICEA import pedigree version")

            # add parameters 1 to 17
            out = ''
            for r in record[0:17]:
                out += r + "\t"

            # Add parameters 18 to 32
            out += \
                brca1t + "\t" + brca1r + "\t" + \
                brca2t + "\t" + brca2r + "\t"

            for r in range(6):  # PALB2, ATM and CHEK2 genetic test results
                out += '0' + "\t"

            out += \
                er + "\t" + \
                pr + "\t" + \
                her2 + "\t" + \
                ck14 + "\t" + \
                ck56 + "\n"
            fout.write(out)

    print("New pedigree file: " + pedigree_file + '.bwa_v4.txt')


class Command(BaseCommand):
    help = 'Convert old style pedigree format to BOADICEA import pedigree v4'

    def add_arguments(self, parser):
        parser.add_argument('file', type=str)

    def handle(self, *args, **options):
        pedigree_file = options['file']
        if not os.path.isfile(pedigree_file):
            raise CommandError('File "%s" does not exist' % pedigree_file)
        if not os.access(pedigree_file, os.R_OK):
            raise CommandError('File "%s" can not be read' % pedigree_file)

        convert_boadicea_pedigree_v4(pedigree_file)
