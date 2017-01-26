from boadicea.exceptions import RiskFactorError
from collections import OrderedDict


class RiskFactors(object):
    ''' Risk factors are defined in terms of a categories. An unknown risk factory
    is category 0, otherwise they are given a non-zero number depending on which
    group they belong to. '''

    categories = OrderedDict([
        ('menarche_age', 7),                # <11, 11, 12, 13, 14, 15, >15
        ('parity', 4),                      # Nulliparous, 1 birth, 2 births, >2 births
        ('age_of_first_live_birth', 4),     # <20, 20-24, 25-29, >29
        ('oral_contraception', 3),          # never, former, current
        ('mht', 4),                         # menopause hormone replacement;
                                            # never, former, current e-type, current c-type
        ('bmi', 4),                         # <18.5, 18.5-24.9, 25-29.9, >=30
        ('alcohol_intake', 7),              # 0g, <5g, 5-14g, 15-24g, 25-34g, 35-44g, >=45g
        ('age_of_menopause', 5),            # <40, 40-44, 45-49, 50-54, >54
        ('mammographic_density', 4)         # Breast Imaging Reporting and Data System;
                                            # BI-RADS 1, 2, 3, 4
        ])

    @staticmethod
    def encode(risk_categories):
        ''' Encode the risk categories into a risk factor. '''
        # Define the number of categories for each factor
        n_categories = list(RiskFactors.categories.values())
        n_factors = len(n_categories)

        # Check that the correct number of command line arguments have been supplied.
        if len(risk_categories) != len(n_categories):
            raise RiskFactorError("Error: Incorrect number of command line arguments specified.\n" +
                                  "This program takes {} arguments, {} supplied".format(len(n_categories),
                                                                                        len(risk_categories)))
        multiplicand = 1
        factor = 0
        for i in range(n_factors):
            # Read in the category for each factor
            try:
                category = int(float(risk_categories[i]))
            except:
                raise RiskFactorError("Error: Unable to convert command line argument number {}, '{}'," +
                                      " to integer.\nThis program takes a list of integers" +
                                      " as arguments".format(i + 1, risk_categories[i]))
            # Check that the category is in bounds.
            if category < 0 or category > n_categories[i]:
                raise RiskFactorError("Error: factor {} out of range, {} > {}".format(i, category, n_categories[i]))

            # Encode the categories into a single factor
            factor += multiplicand * category
            multiplicand = multiplicand * (n_categories[i] + 1)
        return factor

    @staticmethod
    def decode(factor):
        ''' Decode the risk factor into the risk categories. '''
        # Define the number of categories for each factor
        n_categories = list(RiskFactors.categories.values())
        n_factors = len(n_categories)

        # Calcaulate the maximum allowed risk factor code.
        max_factor = 1
        for i in n_categories:
            max_factor *= i + 1
        max_factor -= 1

        # Read in the risk factor code and convert it to integer
        if not isinstance(factor, int):
            raise RiskFactorError("Error: Unable to convert command line argument, {} to integer.\n" +
                                  "This program takes a single integer as argument".format(factor))

        # Check that the category is in bounds
        if factor < 0:
            raise RiskFactorError("Error: factor out of range, {} < {}".format(factor, 0))
        elif factor > max_factor:
            raise RiskFactorError("Error: factor out of range, {} > {}".format(factor, max_factor))

        # Decode the single factor
        dividend = factor
        category = []
        for i in range(n_factors):
            category.append(int(dividend % (n_categories[i] + 1)))
            dividend = (dividend - category[-1]) / (n_categories[i] + 1)
        return category
