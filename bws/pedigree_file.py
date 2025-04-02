"""
Pedigree data file

© 2023 University of Cambridge
SPDX-FileCopyrightText: 2023 University of Cambridge
SPDX-License-Identifier: GPL-3.0-or-later
"""

import logging

from django.conf import settings
from django.utils.translation import gettext_lazy as _

import bws.consts as consts
from bws.exceptions import PedigreeFileError
from bws.pedigree import BwaPedigree, CanRiskPedigree, Pedigree
from bws.risk_factors.bc import BCRiskFactors
from bws.risk_factors.mdensity import Birads, Volpara, Stratus
from bws.risk_factors.oc import OCRiskFactors
from bws.risk_factors.ethnicity import ONSEthnicity
import re


logger = logging.getLogger(__name__)


class Prs(object):
    ''' Polygenic risk score - alpha and zscore values. '''
    def __init__(self, alpha, zscore):
        self.alpha = alpha
        self.zscore = zscore


class CanRiskHeader():
    '''
    CanRisk File Format Header
    '''
    def __init__(self):
        self.lines = []

    def add_line(self, line):
        ''' Append header line to the list of lines. '''
        self.lines.append(line)

    def get_prs(self, val):
        ''' Get a Prs object from a value containing e.g. alpha=float, zscore=float'''
        alpha = zscore = None
        parts = val.replace(" ", "").split(",")
        for p1 in parts:
            p2 = p1.split('=')
            if len(p2) == 2:
                if 'alpha' in p2[0]:
                    alpha = float(p2[1])
                if 'zscore' in p2[0]:
                    zscore = float(p2[1])
                if 'beta' in p2[0]:       # deprecated
                    zscore = float(p2[1])
        if alpha is not None and zscore is not None:
            return Prs(alpha, zscore)
        return None

    def get_risk_factor_codes(self):
        ''' Get breast and ovarian cancer risk factor code, height, BIRADS and PRS from header lines. '''
        bc_rfs = BCRiskFactors()
        oc_rfs = OCRiskFactors()
        bc_prs = oc_prs = pc_prs = None
        hgt = -1
        ons_ethnicity = None
        biobank_ethnicity = None
        md = None
        menopause_status = "0"
        for line in self.lines:
            try:
                parts = line.split('=', 1)
                rfnam = parts[0][2:].lower().strip()    # risk factor name
                rfval = parts[1].strip()                # risk factor value
                if rfnam == 'prs_oc':                   # get ovarian cancer prs
                    oc_prs = self.get_prs(rfval)
                elif rfnam == 'prs_bc':                 # get breast cancer prs
                    bc_prs = self.get_prs(rfval)
                elif rfnam == 'prs_pc':                 # get prostate cancer prs
                    pc_prs = self.get_prs(rfval)
                else:                                   # lookup breast/ovarian cancer risk factors
                    if rfnam == 'height':
                        if rfval == 'NA':
                            continue
                        hgt = float(rfval)
                    elif rfnam == 'birads':
                        md = Birads(rfval)
                    elif rfnam == 'stratus':
                        md = Stratus(rfval)
                    elif rfnam == 'volpara':
                        md = Volpara(rfval)
                    elif rfnam == 'ethnicity':
                        e = rfval.split(';')
                        ons_ethnicity = ONSEthnicity(e[0], e[1] if len(e) > 1 and e[1] != "" else None)
                        biobank_ethnicity = ONSEthnicity.ons2UKBioBank(ons_ethnicity)
                    elif rfnam == 'menopause':
                        menopause_status = "N" if rfval == "N" else "Y"

                    bc_rfs.add_category(rfnam, rfval)
                    oc_rfs.add_category(rfnam, rfval)
            except Exception as e:
                logger.error("CanRisk header format contains an error.", e)
                raise PedigreeFileError("CanRisk header format contains an error in: "+line)

        # add menopause status to volpara/stratus
        if md is not None  and (isinstance(md, Volpara) or isinstance(md, Stratus)):
            md.set_menopause_status(menopause_status)
        return (BCRiskFactors.encode(bc_rfs.cats), OCRiskFactors.encode(oc_rfs.cats), hgt, md, ons_ethnicity, biobank_ethnicity, bc_prs, oc_prs, pc_prs)


class PedigreeFile(object):
    """
    CanRisk and BOADICEA import pedigree file.
    """
    def __init__(self, pedigree_data):
        self.pedigree_data = pedigree_data
        pedigrees_records = [[]]
        canrisk_headers = []
        canrisk_header = CanRiskHeader()
        pid = 0
        famid = None
        file_type = None
        bc_rfc = oc_rfc = 0
        bc_prs = oc_prs = None

        for idx, line in enumerate(pedigree_data.splitlines()):
            if idx == 0:
                if consts.REGEX_CANRISK1_PEDIGREE_FILE_HEADER.match(line):
                    file_type = 'canrisk1'
                    nfields = settings.BOADICEA_CANRISK_FORMAT_ONE_DATA_FIELDS
                elif consts.REGEX_CANRISK2_PEDIGREE_FILE_HEADER.match(line):
                    file_type = 'canrisk2'
                    nfields = settings.BOADICEA_CANRISK_FORMAT_TWO_DATA_FIELDS
                elif consts.REGEX_CANRISK3_PEDIGREE_FILE_HEADER.match(line):
                    file_type = 'canrisk3'
                    nfields = settings.BOADICEA_CANRISK_FORMAT_TWO_DATA_FIELDS
                elif consts.REGEX_CANRISK4_PEDIGREE_FILE_HEADER.match(line):
                    file_type = 'canrisk4'
                    nfields = settings.BOADICEA_CANRISK_FORMAT_FOUR_DATA_FIELDS
                elif consts.REGEX_BWA_PEDIGREE_FILE_HEADER_ONE.match(line):
                    file_type = 'bwa'
                    nfields = settings.BOADICEA_PEDIGREE_FORMAT_FOUR_DATA_FIELDS
                else:
                    raise PedigreeFileError(
                        "The first header record in the pedigree file has unexpected characters. " +
                        "The first header record must be '##CanRisk 3.0'.")
            elif (idx == 1 and file_type == 'bwa') or line.startswith('##FamID'):
                self.column_names = line.replace("##FamID", "FamID").split()
                if (((self.column_names[0] != 'FamID') or
                     (self.column_names[2] != 'Target') or
                     (self.column_names[3] != 'IndivID') or
                     (self.column_names[4] != 'FathID') or
                     (self.column_names[5] != 'MothID'))):
                    raise PedigreeFileError(
                        "Column headers in the pedigree file contains unexpected characters. " +
                        "It must include the 'FamID', 'Target', 'IndivID','FathID' and 'MothID' " +
                        "in columns 1, 3, 4, 5 and 6 respectively.")
            elif line.startswith('##'):
                if '=' in line:                     # risk factor declaration line
                    canrisk_header.add_line(line)
            elif consts.BLANK_LINE.match(line):
                continue
            else:
                delim = ("\t" if line.count("\\t") == nfields-1 else r'\s+')
                record = re.split(delim, line.rstrip())
                if famid is None or famid != record[0]:         # start of pedigree
                    canrisk_headers.append(canrisk_header)
                    canrisk_header = CanRiskHeader()
                if famid is not None and famid != record[0]:    # start of next pedigree found
                    pedigrees_records.append([])
                    pid += 1
                famid = record[0]

                if file_type == 'bwa' and len(record) != nfields:
                    raise PedigreeFileError("A data record has an unexpected number of data items. " +
                                            "BOADICEA format 4 pedigree files should have " +
                                            str(nfields) +
                                            " data items per line.")
                elif file_type == 'canrisk1' and len(record) != nfields:
                    raise PedigreeFileError("A data record has an unexpected number of data items. " +
                                            "CanRisk format 1 pedigree files should have " +
                                            str(nfields) +
                                            " data items per line.")
                elif ((file_type == 'canrisk2' or file_type == 'canrisk3') and len(record) != nfields):
                    raise PedigreeFileError("A data record has an unexpected number of data items. " +
                                            "CanRisk format 2 pedigree files should have " +
                                            str(settings.BOADICEA_CANRISK_FORMAT_TWO_DATA_FIELDS) +
                                            " data items per line.")
                elif ((file_type == 'canrisk4') and len(record) != nfields):
                    raise PedigreeFileError("A data record has an unexpected number of data items. " +
                                            "CanRisk format 4 pedigree files should have " +
                                            str(nfields) +
                                            " data items per line.")
                pedigrees_records[pid].append(line)

        self.pedigrees = []
        for i in range(pid+1):
            if file_type == 'bwa':
                self.pedigrees.append(BwaPedigree(pedigree_records=pedigrees_records[i], file_type=file_type))
            elif file_type.startswith('canrisk'):
                bc_rfc, oc_rfc, hgt, mdensity, ons_ethnicity, biobank_ethnicity, bc_prs, oc_prs, pc_prs = canrisk_headers[i].get_risk_factor_codes()
                self.pedigrees.append(
                    CanRiskPedigree(pedigree_records=pedigrees_records[i], file_type=file_type, delim=delim,
                                    bc_risk_factor_code=bc_rfc, oc_risk_factor_code=oc_rfc,
                                    bc_prs=bc_prs, oc_prs=oc_prs, pc_prs=pc_prs,
                                    hgt=hgt, mdensity=mdensity, ons_ethnicity=ons_ethnicity,
                                    biobank_ethnicity=biobank_ethnicity))

    @classmethod
    def get_incomplete_age_yob(cls, pedigrees):
        if isinstance(pedigrees, Pedigree):
            pedigrees = [pedigrees]
        incomplete = []
        for pedigree in pedigrees:
            people = pedigree.people
            for p in people:
                if not p.is_complete():
                    incomplete.append(p.pid)

        warnings = []
        if len(incomplete) > 0:
            warnings.append(_('year of birth and age at last follow up must be specified in order for ' +
                                      '%(id)s to be included in a calculation') % {'id': ', '.join(incomplete)})
        return warnings
