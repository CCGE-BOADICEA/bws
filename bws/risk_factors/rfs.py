import re
from bws.exceptions import RiskFactorError
import operator


class RiskFactor(object):

    @classmethod
    def camel_to_space(cls, label):
        ''' Convert camel case to space separated. '''
        return re.sub(r'((?<=[a-z])[A-Z]|(?<!\A)[A-Z](?=[a-z]))', r' \1', label)

    @classmethod
    def camel_to_snake(cls, label):
        ''' Convert camel to snake case e.g. AgeOfFirstLiveBirth -> age_of_first_live_birth '''
        return re.sub(r'((?<=[a-z])[A-Z]|(?<!\A)[A-Z](?=[a-z]))', r'_\1', label).lower()

    @classmethod
    def snake_name(cls):
        ''' Return the class name in snake format (e.g. menarche_age) '''
        return cls.camel_to_snake(cls.__name__)

    @classmethod
    def space_name(cls):
        ''' Return the class name in space separated format (e.g. Menarche Age) '''
        return cls.camel_to_space(cls.__name__)

    OPS = {
        ">=": operator.ge,
        ">": operator.gt,
        "<=": operator.le,
        "<": operator.lt,
        "=": operator.eq
    }

    @classmethod
    def eval(cls, val, expr, isreal):
        ''' Evaluate a value against an operand at the start of a string, e.g. <3 '''
        if expr.startswith('<=') or expr.startswith('>='):
            return cls.OPS[expr[0:2]](val, cls.get_num(expr[2:], isreal))
        elif expr.startswith('<') or expr.startswith('>'):
            return cls.OPS[expr[0:1]](val, cls.get_num(expr[1:], isreal))
        else:
            try:
                return cls.OPS["="](val, cls.get_num(expr, isreal))
            except ValueError:
                return False

    @classmethod
    def get_category(cls, val, isreal=False):
        ''' Get category for risk factor. This assumes the categories a'''
        if isinstance(val, str):
            val = val.strip()
        if val == 'NA' or val == '-':
            return 0
        try:
            val = cls.get_num(val, isreal)
            for idx, cat in enumerate(cls.cats):
                if cat == '-':
                    continue

                if '-' in cat:
                    rng = cat.split("-")
                    if rng[0][0] not in cls.OPS:
                        rng[0] = ">="+rng[0]
                    if rng[1][0] not in cls.OPS:
                        rng[1] = "<="+rng[1]
                    if cls.eval(val, rng[0], isreal) and cls.eval(val, rng[1], isreal):
                        return idx
                elif cls.eval(val, cat, isreal):
                    return idx
        except Exception as e:
            print(e)
            raise RiskFactorError("Unknown category for: "+cls.__name__)

    @classmethod
    def get_num(cls, val, isreal):
        if isreal:
            return float(val)
        else:
            return int(val)

    @classmethod
    def isclass(cls, rfname):
        ''' Given a risk factor name determine if it matches the class names or synonym. '''
        return (rfname == cls.__name__.lower() or
                rfname == cls.snake_name().lower() or
                (hasattr(cls, 'synonyms') and rfname in cls.synonyms))

    @classmethod
    def get_value(cls, cat):
        ''' Given a category, e.g. <5 or 14.2-15.4, return a valid value within the category '''
        if cat == '-':
            return cat
        if '-' in cat:
            rng = cat.split("-")
            try:
                float(rng[0])
                return rng[0]
            except Exception:
                return cat

        if cat[0] in cls.OPS:
            if cat[0:2] in cls.OPS:
                return cat[2:]
            try:
                if cat[0] == '<':
                    return str(int(cat[1:])-1)
                elif cat[0] == '>':
                    return str(int(cat[1:])+1)
            except Exception:
                return cat
        return cat


class RiskFactors(object):
    ''' Each risk factor for an individual is defined in terms of a category they are in.
        If a factor is unobserved, missing or not applicable, it is assigned category 0,
        and is not taken into account in the calculation. Otherwise a non-zero number is given
        depending on which group they belong to. These are then combined into a single
        risk factor code (see encode() function) that is used by the BOADICEA risk calculation. '''

    def __init__(self):
        self.cats = [0 for _k in self.categories.keys()]

    @classmethod
    def encode(cls, risk_categories):
        ''' Encode the risk categories into a risk factor. '''
        # Define the number of categories for each factor
        n_categories = list(cls.categories.values())
        n_factors = len(n_categories)

        # Check that the correct number of command line arguments have been supplied.
        if len(risk_categories) != len(n_categories):
            raise RiskFactorError("Incorrect number of risk factors specified.\n" +
                                  "Expecting {} risk factors, {} supplied.".format(len(n_categories),
                                                                                   len(risk_categories)))
        multiplicand = 1
        factor = 0
        for i in range(n_factors):
            # Read in the category for each factor
            try:
                category = int(float(risk_categories[i]))
            except Exception:
                raise RiskFactorError("The {} category '{}' cannot be converted to an integer.".format(
                    cls.risk_factors[i].space_name(), risk_categories[i]))
            # Check that the category is in bounds.
            if category < 0 or category > n_categories[i]:
                raise RiskFactorError("Risk factor ({}) out of range, {} > {}".format(cls.risk_factors[i].space_name(),
                                                                                      category, n_categories[i]))

            # Encode the categories into a single factor
            factor += multiplicand * category
            multiplicand = multiplicand * (n_categories[i] + 1)
        return factor

    @classmethod
    def decode(cls, factor):
        ''' Decode the risk factor into the risk categories. '''
        # Define the number of categories for each factor
        n_categories = list(cls.categories.values())
        n_factors = len(n_categories)
        max_factor = cls.get_max_factor()

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

    @classmethod
    def get_max_factor(cls):
        ''' Calcaulate the maximum allowed risk factor code. '''
        n_categories = list(cls.categories.values())
        max_factor = 1
        for i in n_categories:
            max_factor *= i + 1
        max_factor -= 1
        return max_factor

    def add_category(self, name, val):
        '''
        Given a risk factor name and value add to the category
        '''
        for idx, rf in enumerate(self.risk_factors):
            if rf.isclass(name):
                self.cats[idx] = rf.get_category(val)
