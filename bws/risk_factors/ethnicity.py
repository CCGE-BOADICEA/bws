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
        super().__init__(ethnicity)

    def validate(self):
        if self.ethnicity not in UKBioBankEthnicty.GROUPS:
            raise Exception(self.ethnicity+" not a UK BioBank ethnic group")

    def get_filename(self):
        return "coeffs-BC_" + UKBioBankEthnicty.GROUPS[self.ethnicity] + ".nml"
