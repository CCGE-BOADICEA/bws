from boadicea.exceptions import RiskFactorError
from collections import OrderedDict
from rest_framework.views import APIView
from rest_framework_xml.renderers import XMLRenderer
from rest_framework.renderers import BrowsableAPIRenderer, JSONRenderer
from bws.throttles import BurstRateThrottle, EndUserIDRateThrottle,\
    SustainedRateThrottle
from rest_framework import serializers, status, permissions
from rest_framework.response import Response
from rest_framework.permissions import IsAuthenticated
from rest_framework.authentication import BasicAuthentication,\
    SessionAuthentication, TokenAuthentication
import re


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


class MenarcheAge(RiskFactor):
    cats = ['-', '<11', '11', '12', '13', '14', '15', '>15']
    help_text = 'Age at First Occurrence of Menstruation'


class Parity(RiskFactor):
    cats = ['-', '0', '1', '2', '>2']
    help_text = 'Number of Births (live and still births)'


class AgeOfFirstLiveBirth(RiskFactor):
    cats = ['-', '<20', '20-24', '25-29', '>29']
    help_text = 'Age of First Live Birth'


class OralContraception(RiskFactor):
    cats = ['-', 'never', 'former', 'current']
    help_text = 'Oral Contraception Usage'


class MHT(RiskFactor):
    ''' Menopause hormone replacement '''
    cats = ['-', 'never/former', 'current e-type', 'current other/unknown type (including combined type)']
    help_text = 'Hormone Replacement Therapy'


class BMI(RiskFactor):
    cats = ['-', '<18.5', '18.5-24.9', '25-29.9', '>=30']
    help_text = 'Body Mass Index'


class AlcoholIntake(RiskFactor):
    cats = ['-', '0', '<5', '5-14', '15-24', '25-34', '35-44', '>=45']
    help_text = 'Alcohol Intake (grams/day)'


class AgeOfMenopause(RiskFactor):
    cats = ['-', '<40', '40-44', '45-49', '50-54', '>54']
    help_text = 'Age of Menopause'


class MammographicDensity(RiskFactor):
    cats = ['-', 'BI-RADS 1', 'BI-RADS 2', 'BI-RADS 3', 'BI-RADS 4']
    help_text = 'Mammographic Denisity'


class Height(RiskFactor):
    cats = ['-', '<150.17', '150.17-158.26', '158.26-165.82', '165.82-173.91', '>173.91']
    help_text = 'Height (cm)'


class RiskFactors(object):
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
        max_factor = RiskFactors.get_max_factor()

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

    @staticmethod
    def get_max_factor():
        ''' Calcaulate the maximum allowed risk factor code. '''
        n_categories = list(RiskFactors.categories.values())
        max_factor = 1
        for i in n_categories:
            max_factor *= i + 1
        max_factor -= 1
        return max_factor


class RiskFactorsInputSerializer(serializers.Serializer):
    ''' Risk Factor Categories '''
    for key, value in RiskFactors.categories.items():
        exec(key.lower() + " = serializers.IntegerField(required=False, "
             "max_value="+str(value)+", min_value=0, default=0)")


class RiskFactorsOutputSerializer(serializers.Serializer):
    """ Boadicea result. """
    factor = serializers.IntegerField(max_value=RiskFactors.get_max_factor())


class CanRiskPermission(permissions.BasePermission):
    message = 'Cancer risk factor permission not granted'

    def has_permission(self, request, view):
        return request.user.has_perm('boadicea_auth.can_risk')


class RiskFactorsView(APIView):
    renderer_classes = (XMLRenderer, JSONRenderer, BrowsableAPIRenderer, )
    serializer_class = RiskFactorsInputSerializer
    authentication_classes = (SessionAuthentication, BasicAuthentication, TokenAuthentication, )
    permission_classes = (IsAuthenticated, CanRiskPermission)
    throttle_classes = (BurstRateThrottle, SustainedRateThrottle, EndUserIDRateThrottle)

    def post(self, request):
        """
        Each risk factor for an individual is defined in terms of a category they are in.
        If a factor is unobserved, missing or not applicable, it is assigned category 0,
        and is not taken into account in the calculation. Otherwise a non-zero number is given
        depending on which group they belong to. These are then combined into a single
        risk factor code that is used by the BOADICEA risk calculation.
        ---
        parameters_strategy: merge
        response_serializer: RiskFactorsOutputSerializer
        parameters:
           - name: menarche_age
             description: age of menarche, categories are <11, 11, 12, 13, 14, 15, >15
             type: integer
           - name: parity
             description: categories are Nulliparous, 1 birth, 2 births, >2 births
             type: integer
           - name: age_of_first_live_birth
             description: categories are <20, 20-24, 25-29, >29
             type: integer
           - name: oral_contraception
             description: categories are never, former, current
             type: integer
           - name: mht
             description: categories are never/former, current e-type, current c-type
             type: integer
           - name: bmi
             description: categories are <18.5, 18.5-24.9, 25-29.9, >=30
             type: integer
           - name: alcohol_intake
             description: categories are 0g, <5g, 5-14g, 15-24g, 25-34g, 35-44g, >=45g
             type: integer
           - name: age_of_menopause
             description: categories are <40, 40-44, 45-49, 50-54, >54
             type: integer
           - name: mammographic_density
             description: categories are BI-RADS 1, 2, 3, 4
             type: integer
           - name: height
             description: categories are <150.17, 150.17-158.26, 158.26-165.82, 165.82-173.91, >173.91
             type: integer

        responseMessages:
           - code: 401
             message: Not authenticated

        consumes:
           - application/json
           - application/xml
        produces: ['application/json', 'application/xml']
        """
        serializer = self.serializer_class(data=request.data)
        if serializer.is_valid(raise_exception=True):
            validated_data = serializer.validated_data
            category = []
            for key in RiskFactors.categories.keys():
                category.append(validated_data.get(key))
            factor_serializer = RiskFactorsOutputSerializer({'factor': RiskFactors.encode(category)})
            return Response(factor_serializer.data)
        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)
