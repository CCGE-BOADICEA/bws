'''
Ethnicity

© 2022 Cambridge University
SPDX-FileCopyrightText: 2022 Cambridge University
SPDX-License-Identifier: GPL-3.0-or-later
'''
from abc import ABC, abstractmethod


class Ethnicity(ABC):
    def __init__(self, ethnicity):
        self.ethnicity = ethnicity.lower()
        self.validate()

    @abstractmethod
    def validate(self): raise NotImplementedError

    @abstractmethod
    def get_filename(self): raise NotImplementedError


class UKBioBankEthnicty(Ethnicity):
    '''
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

    @classmethod
    def factory(cls, ethnicity):
        '''
        Parse the ethnicity value to create Ethnicity object
        Examples:
        ##Ethnicity=White;English/Welsh/Scottish/Northern Irish/British;
        ##Ethnicity=White;Irish;
        ##Ethnicity=Mixed/Multiple ethnic groups;White and Black African;
        ##Ethnicity=Asian or Asian British;Chinese;
        ##Ethnicity=Black or Black British;African;
        ##Ethnicity=Other ethnic group;Arab;
        ##Ethnicity=Unknown;
        '''
        parts = ethnicity.split(';')
        e = parts[0].lower().strip()
            
        if e == "white":
            return UKBioBankEthnicty(e);
        elif e == "mixed/multiple ethnic groups": 
            return UKBioBankEthnicty("mixed");#
        elif e == "asian or asian british":
            if parts[1].lower().strip() == "chinese":
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
        
        
    def get_filename(self):
        return UKBioBankEthnicty.GROUPS[self.ethnicity] + ".nml"
