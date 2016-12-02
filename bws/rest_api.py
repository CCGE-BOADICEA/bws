''' API for the BWS REST resources. '''
import datetime
import logging
import shutil
import tempfile

from django.conf import settings
from rest_framework import serializers, status
from rest_framework.authentication import BasicAuthentication, \
    TokenAuthentication, SessionAuthentication

from rest_framework.permissions import IsAuthenticated
from rest_framework.renderers import BrowsableAPIRenderer, JSONRenderer
from rest_framework.response import Response
from rest_framework.throttling import UserRateThrottle, SimpleRateThrottle
from rest_framework.views import APIView
from rest_framework_xml.renderers import XMLRenderer
from boadicea.pedigree import PedigreeFile
from boadicea.calcs import Predictions
from rest_framework.exceptions import NotAcceptable
from rest_framework.compat import is_authenticated
from django.core.files.base import File
# from boadicea.decorator import profile


logger = logging.getLogger(__name__)


class LogThrottleMixin():
    """ Add logging for occasions when throttle stops request. """

    def allow_request(self, request, view):
        """ Implement the check to see if the request should be throttled. """
        success = super(LogThrottleMixin, self).allow_request(request, view)
        if not success:
            return self.throttle_fail(request.user)
        return success

    def throttle_fail(self, user):
        """ Called when a request to the API has failed due to throttling. """
        logger.warning("Request throttled (" + self.__class__.__name__ + "); USER: "+str(user))
        return False


class BurstRateThrottle(LogThrottleMixin, UserRateThrottle):
    """ Throttle short burst of requests from a user. """
    scope = 'burst'


class SustainedRateThrottle(LogThrottleMixin, UserRateThrottle):
    """ Throttle sustained requests from a user. """
    scope = 'sustained'


class EndUserIDRateThrottle(LogThrottleMixin, SimpleRateThrottle):
    """
    Limits the rate of API calls that may be made by a given end user.
    The end user id plus the user will be used as a unique cache key.
    """
    scope = 'enduser_burst'

    def get_cache_key(self, request, view):
        if is_authenticated(request.user):
            ident = str(request.user.pk)
        else:
            ident = str(self.get_ident(request))

        try:
            ident = ident+"$"+request.data.get('user_id')
            return self.cache_format % {
                'scope': self.scope,
                'ident': ident
            }
        except TypeError:
            return None


class PedigreeField(serializers.Field):
    """
    Pedigree field object serialized into a string representation. The field can
    be a str or and uploaded file type.
    """
    def to_representation(self, obj):
        return obj

    def to_internal_value(self, obj):
        assert(isinstance(obj, str) or isinstance(obj, File))

        if isinstance(obj, File):
            pedigree_data = ''
            for chunk in obj.chunks():
                pedigree_data += chunk.decode("utf-8")
            return pedigree_data
        else:
            return obj


class BwsInputSerializer(serializers.Serializer):
    ''' Boadicea result. '''
    user_id = serializers.CharField(min_length=4, max_length=40, required=True)
    pedigree_data = PedigreeField()
    mut_freq = serializers.ChoiceField(choices=['UK', 'Ashkenazi', 'Iceland', 'Custom'],
                                       default='UK', help_text="Mutation frequency")

    MIN_MUT_FREQ = str(settings.MIN_MUTATION_FREQ)
    MAX_MUT_FREQ = str(settings.MAX_MUTATION_FREQ)
    for gene in settings.GENES:
        exec(gene.lower() + "_mut_frequency = serializers.FloatField(required=False, "
             "max_value="+MAX_MUT_FREQ+", min_value="+MIN_MUT_FREQ+")")

    for gene in settings.GENES:
        exec(gene.lower() + "_mut_sensitivity = serializers.FloatField(required=False, default=" +
             str(settings.GENETIC_TEST_SENSITIVITY[gene]) + ", max_value=1, min_value=0)")
    cancer_rates = serializers.ChoiceField(choices=list(settings.CANCER_RATES.keys()))

#     def validate(self, attrs):
#         """ Validate input parameters. """
#         mut_freq = attrs.get('mut_freq')
#         errs = []
#         if mut_freq == 'Custom':
#             for gene in settings.GENES:
#                 mf = attrs.get(gene.lower() + '_mut_frequency')
#                 if mf > settings.MAX_MUTATION_FREQ:
#                     errs.append(gene + " mutation frequency should be less than or equal to " +
#                                 str(settings.MAX_MUTATION_FREQ))
#                 elif mf < settings.MIN_MUTATION_FREQ:
#                     errs.append(gene + " mutation frequency should be greater than or equal to " +
#                                 str(settings.MIN_MUTATION_FREQ))
#         if len(errs) > 0:
#             raise serializers.ValidationError(errs)
#         return attrs
#
#     def isfloat(self, value):
#         """
#         Return true if the given value a float.
#         """
#         try:
#             float(value)
#             return True
#         except:
#             return False


class PedigreeResultSerializer(serializers.Serializer):
    family_id = serializers.CharField(read_only=True)
    cancer_risks = serializers.ListField(read_only=True, required=False)
    mutation_probabilties = serializers.ListField(read_only=True)


class BwsOutputSerializer(serializers.Serializer):
    """ Boadicea result. """
    version = serializers.CharField(read_only=True)
    timestamp = serializers.DateTimeField(read_only=True)
    mutation_frequency = serializers.DictField(read_only=True)
    mutation_sensitivity = serializers.DictField(read_only=True)
    cancer_incidence_rates = serializers.CharField(read_only=True)
    pedigree_result = PedigreeResultSerializer(read_only=True, many=True)
    warnings = serializers.ListField(read_only=True, required=False)


# class CsrfExemptSessionAuthentication(SessionAuthentication):
#
#     def enforce_csrf(self, request):
#         return  # To not perform the csrf check previously happening


class BwsView(APIView):
    renderer_classes = (XMLRenderer, JSONRenderer, BrowsableAPIRenderer, )
    serializer_class = BwsInputSerializer
    authentication_classes = (SessionAuthentication, BasicAuthentication, TokenAuthentication, )
    permission_classes = (IsAuthenticated,)
    throttle_classes = (BurstRateThrottle, SustainedRateThrottle, EndUserIDRateThrottle)

    def get_serializer_class(self):
        return BwsInputSerializer

    # @profile("profile_bws.profile")
    def post(self, request):
        """
        BOADICEA Web-Service (BWS) used to calculate the risks of breast and ovarian cancer and the probability
        that an individual is a carrier of cancer-associated mutations in genes (BRCA1, BRCA2, PALB2, CHEK2, ATM).
        As well as the individuals pedigree, the prediction model takes as input mutation frequency and sensitivity
        for each the genes and the population to use for cancer incidence rates.
        ---
        parameters_strategy: merge
        response_serializer: BwsOutputSerializer
        parameters:
           - name: user_id
             description: unique end user ID, e.g. IP address
             type: string
             required: true
           - name: pedigree_data
             description: BOADICEA pedigree data file
             type: file
             required: true
           - name: mut_freq
             description: mutation frequency
             required: true
             type: string
             paramType: form
             defaultValue: 'UK'
             enum: ['UK', 'Ashkenazi', 'Iceland', 'Custom']
           - name: brca1_mut_frequency
             description: BRCA1 mutation frequency (only available with mut_freq=Custom)
             required: false
             type: float
             paramType: form
           - name: brca2_mut_frequency
             description: BRCA2 mutation frequency (only available with mut_freq=Custom)
             required: false
             type: float
             paramType: form
           - name: palb2_mut_frequency
             description: PALB2 mutation frequency (only available with mut_freq=Custom)
             required: false
             type: float
             paramType: form
           - name: atm_mut_frequency
             description: ATM mutation frequency (only available with mut_freq=Custom)
             required: false
             type: float
             paramType: form
           - name: chek2_mut_frequency
             description: CHEK2 mutation frequency (only available with mut_freq=Custom)
             required: false
             type: float
             paramType: form
           - name: cancer_rates
             description: cancer incidence rates
             required: true
             type: string
             paramType: form
             defaultValue: 'UK'
             enum: ['UK', 'UK-version-1', 'Australia', 'USA-white', 'Denmark', 'Finland',
             'Iceland', 'New-Zealand', 'Norway', 'Sweden']

           - name: brca1_mut_sensitivity
             description: BRCA1 mutation sensitivity
             required: false
             type: float
             paramType: form
             defaultValue: 0.7
           - name: brca2_mut_sensitivity
             description: BRCA2 mutation sensitivity
             required: false
             type: float
             paramType: form
             defaultValue: 0.8
           - name: palb2_mut_sensitivity
             description: PALB2 mutation sensitivity
             required: false
             type: float
             paramType: form
             defaultValue: 0.9
           - name: atm_mut_sensitivity
             description: ATM mutation sensitivity
             required: false
             type: float
             paramType: form
             defaultValue: 0.9
           - name: chek2_mut_sensitivity
             description: CHEK2 mutation sensitivity
             required: false
             type: float
             paramType: form
             defaultValue: 1.0

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
            pedigree_data = validated_data.get('pedigree_data')

            pf = PedigreeFile(pedigree_data)
            population = validated_data.get('mut_freq', 'UK')
            cancer_rates = settings.CANCER_RATES.get(validated_data.get('cancer_rates'))

            if population != 'Custom':
                mutation_frequency = settings.MUTATION_FREQUENCIES[population]
            else:
                mutation_frequency = {}
                for gene in settings.GENES:
                    try:
                        mutation_frequency[gene] = float(validated_data.get(gene.lower() + '_mut_frequency'))
                    except TypeError:
                        raise NotAcceptable("Invalid mutation frequency for " + gene + ".")

            mutation_sensitivity = {
                k: float(validated_data.get(k.lower() + "_mut_sensitivity", settings.GENETIC_TEST_SENSITIVITY[k]))
                for k in settings.GENETIC_TEST_SENSITIVITY.keys()
            }

            output = {
                "version": settings.BOADICEA_VERSION,
                "timestamp": datetime.datetime.now(),
                "mutation_frequency": {population: mutation_frequency},
                "mutation_sensitivity": mutation_sensitivity,
                "cancer_incidence_rates": cancer_rates,
                "pedigree_result": []
            }

            warnings = PedigreeFile.validate(pf.pedigrees)
            if len(warnings) > 0:
                output['warnings'] = warnings

            cwd = tempfile.mkdtemp(prefix=str(request.user)+"_", dir="/tmp")
            try:
                for pedi in pf.pedigrees:
                    this_pedigree = {}
                    this_pedigree["family_id"] = pedi.famid

                    calcs = Predictions(pedi, mutation_frequency=mutation_frequency,
                                        mutation_sensitivity=mutation_sensitivity,
                                        cancer_rates=cancer_rates, cwd=cwd, request=request)
                    try:
                        if calcs.mutation_probabilties is not None:
                            this_pedigree["mutation_probabilties"] = calcs.mutation_probabilties
                    except AttributeError as e:
                        if 'warnings' in output:
                            output['warnings'].append('mutation probabilties not provided')
                        else:
                            output['warnings'] = 'mutation probabilties not provided'
                        logger.debug('cancer risks not provided :: '+str(e))
                    try:
                        if calcs.cancer_risks is not None:
                            this_pedigree["cancer_risks"] = calcs.cancer_risks
                    except AttributeError as e:
                        if 'warnings' in output:
                            output['warnings'].append('cancer risks not provided')
                        else:
                            output['warnings'] = 'cancer risks not provided'
                        logger.debug('cancer risks not provided :: '+str(e))

                    output["pedigree_result"].append(this_pedigree)
            finally:
                shutil.rmtree(cwd)
            output_serialiser = BwsOutputSerializer(output)
            return Response(output_serialiser.data)
        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)
