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
from bws.risk_factors.bc import BCRiskFactors


class RiskFactorsInputSerializer(serializers.Serializer):
    ''' Risk Factor Categories '''
    for key, value in BCRiskFactors.categories.items():
        exec(key.lower() + " = serializers.IntegerField(required=False, "
             "max_value="+str(value)+", min_value=0, default=0)")


class RiskFactorsOutputSerializer(serializers.Serializer):
    """ Boadicea result. """
    factor = serializers.IntegerField(max_value=BCRiskFactors.get_max_factor())


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
            for key in BCRiskFactors.categories.keys():
                category.append(validated_data.get(key))
            factor_serializer = RiskFactorsOutputSerializer({'factor': BCRiskFactors.encode(category)})
            return Response(factor_serializer.data)
        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)
