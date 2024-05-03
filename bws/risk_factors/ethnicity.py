'''
Ethnicity

© 2023 University of Cambridge
SPDX-FileCopyrightText: 2023 University of Cambridge
SPDX-License-Identifier: GPL-3.0-or-later
'''
from abc import ABC, abstractmethod


class Ethnicity(ABC):
    def __init__(self, ethnicity, ethnicityBackground=None):
        self.ethnicity = ethnicity.lower()
        self.ethnicityBackground = ethnicityBackground.lower() if ethnicityBackground is not None else None
        self.validate()

    @abstractmethod
    def validate(self): raise NotImplementedError

    @abstractmethod
    def get_filename(self): raise NotImplementedError


class ONSEthnicity(Ethnicity):
    '''
    https://www.ons.gov.uk/methodology/classificationsandstandards/measuringequality/ethnicgroupnationalidentityandreligion#ethnic-group
    '''
    GROUPS = {
        "White" :
            ["English/Welsh/Scottish/Northern Irish/British",
             "Irish",
             "Gypsy or Irish Traveller",
             "Any other White background, please describe"
            ],
        "Mixed/Multiple ethnic groups":
            ["White and Black Caribbean",
             "White and Black African", 
             "White and Asian",
             "Any other Mixed/Multiple ethnic background, please describe"
            ],
        "Asian or Asian British":
            ["Indian",
             "Pakistani",
             "Bangladeshi",
             "Chinese",
             "Any other Asian background, please describe"
            ],
        "Black or Black British":
            ["African",
             "Caribbean",
             "Any other Black/African/Caribbean background, please describe"
            ],
        "Other ethnic group":
            ["Arab",
             "Any other ethnic group, please describe"
            ],
        "Unknown" : None
    }
    
    # GROUPS_LOWERCASE = dict((k.lower(), (v.lower() for v in vs) if vs is not None else None) for k,vs in GROUPS.items())
    GROUPS_LOWERCASE = dict((k.lower(), ([v.lower() for v in vs]) if vs is not None else None) for k,vs in GROUPS.items())

    def __init__(self, ethnicity="Unknown", ethnicityBackground=None):
        assert(ethnicity in ONSEthnicity.GROUPS)
        super().__init__(ethnicity, ethnicityBackground)

    def validate(self):
        if self.ethnicity not in ONSEthnicity.GROUPS_LOWERCASE:
            raise Exception(self.ethnicity.title()+" not an ONS ethnic group")
        if (self.ethnicityBackground is not None and 
            self.ethnicityBackground not in ONSEthnicity.GROUPS_LOWERCASE[self.ethnicity]):
            raise Exception(self.ethnicityBackground+" not an ethnic background for the ONS ethnic group: "+self.ethnicity)

    def get_filename(self): raise NotImplementedError

    ''' Get string representation '''
    def get_string(self):
        if self.ethnicityBackground is not None:
            return self.ethnicity + ";" + self.ethnicityBackground
        return self.ethnicity

    @classmethod
    def ons2UKBioBank(cls, onsEthnicity):
        '''
        Map ONS to UK BioBank ethnicity
        '''
        e = onsEthnicity.ethnicity
        ebg = onsEthnicity.ethnicityBackground
            
        if e == "white":
            return UKBioBankEthnicty(e);
        elif e == "mixed/multiple ethnic groups": 
            return UKBioBankEthnicty("mixed");#
        elif e == "asian or asian british":
            if ebg == "chinese":
                return UKBioBankEthnicty("chinese")
            else:
                return UKBioBankEthnicty("asian")
        elif e ==  "black or black british":
            return UKBioBankEthnicty("black")
        elif e ==  "other ethnic group":
            return UKBioBankEthnicty("other")
        elif e ==  "unknown":
            return UKBioBankEthnicty("unknown")            
        return UKBioBankEthnicty()

class UKBioBankEthnicty(Ethnicity):
    '''
    https://biobank.ndph.ox.ac.uk/ukb/field.cgi?id=21000
    
    UK Biobank ethnic groups:
        1 : White
        2 : Mixed
        3 : Asian or Asian British
        4 : Black or Black British
        5 : Chinese
        6 : Other ethnic group
        7 : Do not know
        8 : Prefer not to answer
    '''

    GROUPS = {
        "white": "UK-european",
        "mixed": "UK-mixed",
        "asian": "UK-southAsian",
        "black": "UK-african",
        "chinese": "UK-eastAsian",
        "other": "UK-other",
        "unknown": "UK-unknown",
        "na": "UK-pop"
    }

    def __init__(self, ethnicity="na"):
        assert(ethnicity in UKBioBankEthnicty.GROUPS)
        super().__init__(ethnicity)

    def validate(self):
        if self.ethnicity not in UKBioBankEthnicty.GROUPS:
            raise Exception(self.ethnicity+" not a UK BioBank ethnic group")        
        
    def get_filename(self):
        return UKBioBankEthnicty.GROUPS[self.ethnicity] + ".nml"

    def get_group(self):
        if self.ethnicity == "na":
            return "UK"
        return "UK "+self.ethnicity.title()
