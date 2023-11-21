'''
Mammographic density - supports BIRADS, Status and Volpara

Encoding for the pedigree (Fortran) file:
If employing BIRAD, the value should be an integer between 1 and 4. If employing
continuous methods, the value should be a real number in the form N.xxxx. 
In this case, N refers to the method (10=Stratus, 20= Volpara) whereas xxxxx
is the mammographic density percentage. 
Example 1: MD = 42.42% measured with Volpara should be coded as “20.42420”
Example 2: MD = category 3 of Birads should be coded as “3”

© 2023 University of Cambridge
SPDX-FileCopyrightText: 2023 University of Cambridge
SPDX-License-Identifier: GPL-3.0-or-later
'''
from abc import ABC, abstractmethod


class MammographicDensity(ABC):
    def __init__(self, md):
        self.md = md
    
    @abstractmethod
    def get_pedigree_str(self): raise NotImplementedError
    
    @abstractmethod
    def get_display_str(self): raise NotImplementedError


class Birads(MammographicDensity):
    
    def get_pedigree_str(self):
        mdcat = Birads.get_category(self.md)
        if mdcat > 0:
            return "0000000"+str(mdcat)
        return "00000000"
    
    def get_display_str(self):
        return "BI-RADS "+str(Birads.get_category(self.md))
        
    @classmethod
    def get_category(cls, val):
        alt1 = ['NA', 'a', 'b', 'c', 'd']
        alt2 = ['NA', '1', '2', '3', '4']
        val = val.lower().replace('bi-rads ', '')
        cats = ['-', 'BI-RADS a', 'BI-RADS b', 'BI-RADS c', 'BI-RADS d']
        for idx, _cat in enumerate(cats):
            if val == alt1[idx] or val == alt2[idx]:
                return idx
        return 0

class Stratus(MammographicDensity):

    def get_pedigree_str(self):
        stratus = 10 + (float(self.md)/100)
        return format(stratus, '.5f')
    
    def get_display_str(self):
        return "Stratus "+str(self.md)


class Volpara(MammographicDensity):

    def get_pedigree_str(self):
        volpara = 20 + (float(self.md)/100)
        return format(volpara, '.5f')
    
    def get_display_str(self):
        return "Volpara "+str(self.md)
