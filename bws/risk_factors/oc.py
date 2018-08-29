'''
Ovarian cancer risk factors.
'''

from bws.risk_factors.rfs import RiskFactor, RiskFactors
from collections import OrderedDict


class Parity(RiskFactor):
    cats = ['-', '0', '1', '>1']
    help_text = 'Number of Births (live and still births)'


class OralContraception(RiskFactor):
    cats = ['-', 'never', '<5', '>=5']
    help_text = 'Oral Contraception Usage'
    synonyms = ['oc_use']

    @classmethod
    def get_category(cls, val, isreal=False):
        if val == '-':
            return 0
        if val == 'N' or val.lower() == 'never':
            return 1
        if ':' not in val:
            return 0
        val = val.split(':')[1]
        return super(OralContraception, cls).get_category(val, isreal)


class MHT(RiskFactor):
    ''' Menopause hormone replacement '''
    cats = ['-', 'never', 'ever']
    help_text = 'Hormone Replacement Therapy'
    synonyms = ['mht_use']

    @classmethod
    def get_category(cls, val):
        if val == 'N':
            return 1
        elif val == 'E' or val == 'F' or val == 'C':
            return 2
        return 0


class TubalLigation(RiskFactor):
    ''' Tubal ligation surgical procedure for sterilization '''
    cats = ['-', 'no', 'yes']
    help_text = 'Tubal ligation surgical procedure for sterilization'
    synonyms = ['TL']

    @classmethod
    def get_category(cls, val):
        alt = ['na', 'n', 'y']
        val = val.lower()
        for idx, cat in enumerate(TubalLigation.cats):
            if val == cat or val == alt[idx]:
                return idx
        return 0


class Endometriosis(RiskFactor):
    cats = []
    synonyms = ['Endo']


class BMI(RiskFactor):
    cats = ['-', '<22.5', '22.5-<25', '25-<27.5', '27.5-<30', '>=30']
    help_text = 'Body Mass Index'

    @classmethod
    def get_category(cls, val, isreal=True):
        return super(BMI, cls).get_category(val, isreal)


class Height(RiskFactor):
    cats = ['-', '<=152.9', '152.91-159.64', '159.65-165.95', '165.96-172.69', '>172.70']
    help_text = 'Height (cm)'

    @classmethod
    def get_category(cls, val, isreal=True):
        return super(Height, cls).get_category(val, isreal)


class OCRiskFactors(RiskFactors):
    ''' Each risk factor for an individual is defined in terms of a category they are in.
        If a factor is unobserved, missing or not applicable, it is assigned category 0,
        and is not taken into account in the calculation. Otherwise a non-zero number is given
        depending on which group they belong to. These are then combined into a single
        risk factor code (see encode() function) that is used by the BOADICEA risk calculation. '''

    # list of risk factor classes
    risk_factors = [
        Parity,
        OralContraception,
        MHT,
        TubalLigation,
        Endometriosis,
        BMI,
        Height
    ]

    # dictionary of risk factor name and number of categories
    categories = OrderedDict(
        (rf.snake_name(), len(rf.cats)-1) for rf in risk_factors
    )

    # dictionary of risk factor name and the categories
    risk_factors_categories = OrderedDict(
        (rf.snake_name(), rf.cats) for rf in risk_factors
    )
