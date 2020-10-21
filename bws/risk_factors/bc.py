'''
Breast cancer risk factors.
'''

from bws.risk_factors.rfs import RiskFactor, RiskFactors
from collections import OrderedDict
from django.utils.translation import gettext_lazy as _


class MenarcheAge(RiskFactor):
    cats = ['-', '<11', '11', '12', '13', '14', '15', '>15']
    help_text = _('Age at First Occurrence of Menstruation')
    synonyms = ['menarche']


class Parity(RiskFactor):
    cats = ['-', '0', '1', '2', '>2']
    help_text = _('Number of Children')


class AgeOfFirstLiveBirth(RiskFactor):
    cats = ['-', '<20', '20-24', '25-29', '>29']
    help_text = _('Age of First Live Birth')
    synonyms = ['first_live_birth']


class OralContraception(RiskFactor):
    cats = ['-', 'never', 'former', 'current']
    help_text = _('Oral Contraception Usage')
    synonyms = ['oc_use']

    @classmethod
    def get_category(cls, val):
        alt = ['na', 'n', 'f', 'c']
        if ':' in val:
            val = val.split(':')[0]
        val = val.lower()
        for idx, cat in enumerate(OralContraception.cats):
            if val == cat or val == alt[idx]:
                return idx
        return 0


class MHT(RiskFactor):
    ''' Menopause hormone replacement '''
    cats = ['-', 'never/former', 'current e-type', 'current other/unknown type (including combined type)']
    help_text = _('Hormone Replacement Therapy')
    synonyms = ['mht_use']

    @classmethod
    def get_category(cls, val):
        if val == 'N' or val == 'F' or val == 'never/former':
            return 1
        elif val == 'E' or val == 'current e-type':
            return 2
        elif val == 'C' or val == 'current other/unknown type (including combined type)':
            return 3
        return 0


class BMI(RiskFactor):
    cats = ['-', '<18.5', '18.5-<25', '25-<30', '>=30']
    help_text = _('Body Mass Index')

    @classmethod
    def get_category(cls, val, isreal=True):
        return super(BMI, cls).get_category(val, isreal)


class AlcoholIntake(RiskFactor):
    cats = ['-', '0', '<5', '5-<15', '15-<25', '25-<35', '35-<45', '>=45']

    help_text = _('Alcohol Intake (grams/day)')
    synonyms = ['alcohol']

    @classmethod
    def get_category(cls, val, isreal=True):
        return super(AlcoholIntake, cls).get_category(val, isreal)


class AgeOfMenopause(RiskFactor):
    cats = ['-', '<40', '40-44', '45-49', '50-54', '>54']
    help_text = _('Age of Menopause')
    synonyms = ['menopause']


class MammographicDensity(RiskFactor):
    cats = ['-', 'BI-RADS a', 'BI-RADS b', 'BI-RADS c', 'BI-RADS d']
    help_text = _('Mammographic Density')
    synonyms = ['birads']

    @classmethod
    def get_category(cls, val):
        alt1 = ['NA', 'a', 'b', 'c', 'd']
        alt2 = ['NA', '1', '2', '3', '4']
        val = val.lower().replace('bi-rads ', '')
        for idx, _cat in enumerate(cls.cats):
            if val == alt1[idx] or val == alt2[idx]:
                return idx
        return 0


class Height(RiskFactor):
    cats = ['-', '<152.91', '152.91-<159.65', '159.65-<165.96', '165.96-<172.70', '>=172.70']
    help_text = _('Height (cm)')

    @classmethod
    def get_category(cls, val, isreal=True):
        return super(Height, cls).get_category(val, isreal)


class BCRiskFactors(RiskFactors):
    ''' Each risk factor for an individual is defined in terms of a category they are in.
        If a factor is unobserved, missing or not applicable, it is assigned category 0,
        and is not taken into account in the calculation. Otherwise a non-zero number is given
        depending on which group they belong to. These are then combined into a single
        risk factor code (see encode() function) that is used by the BOADICEA risk calculation. '''

    # list of risk factor classes
    risk_factors = [
        MenarcheAge,
        Parity,
        AgeOfFirstLiveBirth,
        OralContraception,
        MHT,
        BMI,
        AlcoholIntake,
        AgeOfMenopause,
        MammographicDensity,
        Height
    ]

    # dictionary of risk factor name and number of categories
    categories = OrderedDict(
        (rf.snake_name(), len(rf.cats)-1) for rf in risk_factors
    )

#     categories = OrderedDict([
#         ('menarche_age', 7),                # <11, 11, 12, 13, 14, 15, >15
#         ('parity', 4),                      # Nulliparous, 1 birth, 2 births, >2 births
#         ('age_of_first_live_birth', 4),     # <20, 20-24, 25-29, >29
#         ('oral_contraception', 3),          # never, former, current
#         ('mht', 4),                         # menopause hormone replacement;
#                                             # never, former, current e-type, current c-type
#         ('bmi', 4),                         # <18.5, 18.5-24.9, 25-29.9, >=30
#         ('alcohol_intake', 7),              # 0g, <5g, 5-14g, 15-24g, 25-34g, 35-44g, >=45g
#         ('age_of_menopause', 5),            # <40, 40-44, 45-49, 50-54, >54
#         ('mammographic_density', 4),        # Breast Imaging Reporting and Data System; BI-RADS 1, 2, 3, 4
#         ('height', 5)                       # Height/cm <150.17, 150.17-158.26, 158.26-165.82, 165.82-173.91, >173.91
#         ])

    # dictionary of risk factor name and the categories
    risk_factors_categories = OrderedDict(
        (rf.snake_name(), rf.cats) for rf in risk_factors
    )
