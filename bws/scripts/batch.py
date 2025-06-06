"""
Utilities to assist in running batch processing script and
extract results.

© 2023 University of Cambridge
SPDX-FileCopyrightText: 2023 University of Cambridge
SPDX-License-Identifier: GPL-3.0-or-later
"""

import os
from subprocess import PIPE, Popen
from bws.pedigree_file import PedigreeFile, CanRiskPedigree
import re
import math
from bws.risk_factors.ethnicity import ONSEthnicity, UKBioBankEthnicty


def run_batch(FORTRAN, cwd, csvfile, ofile, irates, ashkn=False, mut_freq="UK", model='BC', muts=False, 
              verbose=False, biobank_ethnicity=UKBioBankEthnicty()):
    ''' Run batch processing script. '''
    if ashkn or mut_freq == "ASHKENAZI":
        setting = FORTRAN+"settings_"+model+"_AJ"+".ini"
    else:
        if model != "BC" or UKBioBankEthnicty.GROUPS[biobank_ethnicity.ethnicity] ==  "UK-pop":
            setting = FORTRAN+"settings_"+model+"_"+mut_freq+".ini"
        else:
            setting = FORTRAN+"settings_"+model+"_UK_non_european.ini"

    if model  == 'OC':
        model_path = os.path.join(FORTRAN, "Ovarian-Cancer-Model")
    elif model == 'BC':
        model_path = os.path.join(FORTRAN, "Breast-Cancer-Model")
    elif model == 'PC':
        model_path = os.path.join(FORTRAN, "Prostate-Cancer-Model")
    cmd = [FORTRAN+"run_job.sh",
           "-r", ofile,
           "-i", os.path.join(model_path, "Data/incidences_"+irates.replace('New-Zealand', 'New_Zealand')),
           "-s", setting,
           "-l", os.path.join(cwd, model+"runlog.log")]
    
    if muts:
        cmd.append('-m')
        cmd.append('p')

    if model == 'OC':
        cmd.append('-c')
        cmd.append('o')
        cmd.append('-e')
        cmd.append(os.path.join(model_path, "Data/coeffs-OC_"+biobank_ethnicity.get_filename()))
    elif model == 'BC':
        cmd.append('-c')
        cmd.append('b')
        cmd.append('-e')
        cmd.append(os.path.join(model_path, "Data/coeffs-BC_"+biobank_ethnicity.get_filename()))
    elif model == 'PC':
        cmd.append('-c')
        cmd.append('p')
        cmd.append('-e')
        cmd.append(os.path.join(model_path, "Data/coeffs-PC_"+biobank_ethnicity.get_filename()))

    cmd.append(csvfile)
    
    if verbose:
        print(' '.join(cmd))

    process = Popen(cmd, cwd=FORTRAN, stdout=PIPE, stderr=PIPE)
    (outs, errs) = process.communicate()
    _exit_code = process.wait()
    if _exit_code != 0:
        print(outs)
        print(errs)
        print(' '.join(cmd))
        print("EXIT CODE "+str(_exit_code))
        exit(1)
    return outs, errs


def get_batch_results(fname, calc_ages):
    ''' Get risk from batch '''
    if not os.path.isfile(fname):
        print("BATCH RESULTS FILE NOT FOUND :: "+fname)
        return None

    c_batch = {}
    f = open(fname, "r")
    for line in f:
        crisks = line.strip().split(',')
        if len(crisks) == 4 and not line.startswith('file') and 'Censor' not in line:
            censoring_age = int(crisks[2])
            if censoring_age in calc_ages:
                c_batch[censoring_age] = crisks[3]
    f.close()
    return c_batch


def get_mp(fname):
    ''' Get mutation propabilities from batch '''
    if not os.path.isfile(fname):
        print("BATCH FILE NOT FOUND :: "+fname)
        return None
    if os.stat(fname).st_size == 0:
        return None

    f = open(fname, "r")
    for line in f:
        if 'PROBABILITIES' in line:
            continue
        parts = line.strip().replace('NO_PATHOGENIC_VARIANTS', 'no mutation').split(',')[2:]
        if 'NO_PATHOGENIC_VARIANTS' in line:
            keys = parts
        else:
            vals = parts
    f.close()
    return {keys[idx]: val for idx, val in enumerate(vals)}


def get_censoring_ages(bwa):
    with open(bwa, 'r') as f:
        pedigree_data = f.read()
    f.close()
    pedigree = PedigreeFile(pedigree_data).pedigrees[0]
    # get ages to calculate risks
    target = pedigree.get_target()
    tage = int(target.age)      # target age at last follow up
    c_ages = []
    alf = tage
    while alf <= 79:
        alf += 1
        if alf-tage in [1, 5, 10]:
            c_ages.append(alf)
    c_ages.append(80)
    return c_ages


def add_prs(line, cancer, rfsnames, rfs):
    '''
    Add PRS arrays for csv batch input file parameters.
    @param line: CanRisk PRS header e.g. PRS_BC=alpha=0.444,zscore=1.12
    @param cancer: string denoting cancer type, i.e. 'BC' or 'OC'
    @param rfsnames: array of risk factor names
    @param rfs: risk factor values
    '''
    zscore = re.match("##PRS.*(zscore=([-]?\d*\.\d+)).*", line)
    alpha = re.match("##PRS.*(alpha=([-]?\d*\.\d+)).*", line)
    if zscore is not None:
        rfsnames.append(['PRS_'+cancer+'_z', cancer+'_PRS_z'])
        rfs['PRS_'+cancer+'_z'] = "%8.5f" % float(zscore.group(2))
    if alpha is not None:
        rfsnames.append(['PRS_'+cancer+'_alpha', cancer+'_PRS_alpha'])
        rfs['PRS_'+cancer+'_alpha'] = alpha.group(2)


def get_menopause_status_code(bwa):
    '''
    Menopause code to be used with Volpara/Stratus mammographic densities
    '''
    f = open(bwa, "r")
    try:
        for line in f:
            if line.startswith("##") and "##CanRisk" not in line and "##FamID" not in line:
                line = line.replace("##", "").strip().split("=")
                name = line[0].capitalize()
                if name == "Menopause":
                    if line[1] == "N":
                        return float(1)
                    else:
                        return float(2)
    finally:
        f.close()
    return float(0)
    
def get_rfs(bwa):
    '''  Get risk factor names and values plus PRS from CanRisk file for CSV file '''
    rfsnames = []
    rfs = {}
    biobank_ethnicity = UKBioBankEthnicty()
    menopause_status_code = get_menopause_status_code(bwa)
    f = open(bwa, "r")    
    for line in f:
        if line.startswith("##") and "##CanRisk" not in line and "##FamID" not in line:
            if "PRS_BC" in line:      # alpha=0.45,zscore=0.1234
                add_prs(line, 'BC', rfsnames, rfs)
            elif "PRS_OC" in line:    # alpha=0.45,zscore=0.1234
                add_prs(line, 'OC', rfsnames, rfs)
            elif "PRS_PC" in line:    # alpha=0.45,zscore=0.1234
                add_prs(line, 'PC', rfsnames, rfs)
            elif "ethnicity" in line.lower():
                line = line.replace("##ethnicity=", "").strip().split(";")
                if len(line)==2 and line[1] == "":
                    line[1] = None
                onsEthnicity = ONSEthnicity(line[0], line[1] if len(line)==2 else None)
                biobank_ethnicity = ONSEthnicity.ons2UKBioBank(onsEthnicity)
            else:
                line = line.replace("##", "").strip().split("=")

                if line[0].isupper():
                    if line[0] == 'TL':
                        name = 'Tubal_Ligation'
                    else:
                        name = line[0]
                elif line[0] == 'mht_use':
                    name = 'MHT_use'
                elif line[0] == 'height':
                    name = line[0]
                    line[1] = ("%8.4f" % float(line[1]))
                elif line[0] == 'oc_use':
                    name = 'OC_Use'
                    if line[1] == "N" or line[1] == "C":
                        rfsnames.append(['OC_Duration', 'OC_Duration'])
                        rfs['OC_Duration'] = '0'
                    elif ":" in line[1]:
                        rfsnames.append(['OC_Duration', 'OC_Duration'])
                        parts = line[1].split(":")
                        rfs['OC_Duration'] = parts[1]
                        line[1] = parts[0]
                elif line[0] == 'birads':
                    name = 'Mamm_density'
                elif line[0] == 'stratus':
                    name = 'Mamm_density'
                    line[1] = str(10+(float(line[1])/100.)+menopause_status_code)
                elif line[0] == 'volpara':
                    name = 'Mamm_density'
                    line[1] = str(20+(float(line[1])/100.)+menopause_status_code)
                elif line[0] == 'endo':
                    name = 'Endometriosis'
                else:
                    name = line[0].capitalize()

                if name == "Menopause" and line[1] == "N":
                    continue            
                rfsnames.append([line[0], name])
                rfs[line[0]] = line[1]

    f.close()

    with open(bwa, 'r') as f:
        pedigree_data = f.read()

    pedigree = PedigreeFile(pedigree_data).pedigrees[0]
    ashkn = pedigree.is_ashkn()
    return rfsnames, rfs, ashkn, biobank_ethnicity


def compare_mp(model, mp_batch, mp_ws, exact_matches, abs_tol=1e-09):
    ''' '''
    exact = True
    msg = ""
    for k, v in mp_batch.items():
        if math.isclose(float(v), float(mp_ws[k]), abs_tol=abs_tol) or k == "no mutation":
            msg += k+":"+mp_ws[k]+"="+v+" "
        else:
            msg += k+":"+mp_ws[k]+"?"+v+" "
            exact_matches += 1
            exact = False
            print(k+": DIFFERENCE ["+str(float(v)-float(mp_ws[k]))+"]")
    if exact:
        print(model+" EXACT MATCH ::: "+msg)
    else:
        print(model+" DIFFERENCE  ::: "+msg)
    return exact_matches


def is_canrisk(bwa):
    '''  Return true if CanRisk file type '''
    with open(bwa, 'r') as f:
        pedigree_data = f.read()

    pedigree = PedigreeFile(pedigree_data).pedigrees[0]
    return isinstance(pedigree, CanRiskPedigree)
