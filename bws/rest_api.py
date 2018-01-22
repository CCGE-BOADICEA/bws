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
from rest_framework.renderers import JSONRenderer, TemplateHTMLRenderer   # , BrowsableAPIRenderer
from rest_framework.response import Response
from rest_framework.views import APIView
from rest_framework_xml.renderers import XMLRenderer
from boadicea.pedigree import PedigreeFile
from boadicea.calcs import Predictions
from rest_framework.exceptions import NotAcceptable, ValidationError
from django.core.files.base import File
from bws.throttles import BurstRateThrottle, EndUserIDRateThrottle,\
    SustainedRateThrottle
from bws.risk_factors import RiskFactors
from django.http.response import JsonResponse

# from boadicea.decorator import profile


logger = logging.getLogger(__name__)


class FileField(serializers.Field):
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
    pedigree_data = FileField()
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


class Prs(object):

    def __init__(self, alpha, beta):
        ''' Polygenic risk alpha and beta values calculated from VCF file. '''
        self.alpha = alpha
        self.beta = beta


class BwsExtendedInputSerializer(BwsInputSerializer):
    ''' Other input parameters. '''
    risk_factor_code = serializers.IntegerField(max_value=RiskFactors.get_max_factor(),
                                                min_value=0, default=0)
    prs = serializers.JSONField(required=False)

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
    proband_id = serializers.CharField(read_only=True)
    cancer_risks = serializers.ListField(read_only=True, required=False)
    baseline_cancer_risks = serializers.ListField(read_only=True, required=False)
    lifetime_cancer_risk = serializers.ListField(read_only=True, required=False)
    baseline_lifetime_cancer_risk = serializers.ListField(read_only=True, required=False)
    ten_yr_cancer_risk = serializers.ListField(read_only=True, required=False)
    baseline_ten_yr_cancer_risk = serializers.ListField(read_only=True, required=False)
    mutation_probabilties = serializers.ListField(read_only=True)


class BwsOutputSerializer(serializers.Serializer):
    """ Boadicea result. """
    version = serializers.CharField(read_only=True)
    timestamp = serializers.DateTimeField(read_only=True)
    mutation_frequency = serializers.DictField(read_only=True)
    mutation_sensitivity = serializers.DictField(read_only=True)
    cancer_incidence_rates = serializers.CharField(read_only=True)
    prs = serializers.DictField(read_only=True, required=False)
    risk_factors = serializers.DictField(read_only=True, required=False)
    pedigree_result = PedigreeResultSerializer(read_only=True, many=True)
    warnings = serializers.ListField(read_only=True, required=False)


class CsrfExemptSessionAuthentication(SessionAuthentication):

    def enforce_csrf(self, request):
        return  # To not perform the csrf check previously happening


class BwsView(APIView):
    renderer_classes = (XMLRenderer, JSONRenderer, TemplateHTMLRenderer, )
    serializer_class = BwsExtendedInputSerializer
    authentication_classes = (CsrfExemptSessionAuthentication, BasicAuthentication, TokenAuthentication, )
    permission_classes = (IsAuthenticated,)
    throttle_classes = (BurstRateThrottle, SustainedRateThrottle, EndUserIDRateThrottle)

    def get_serializer_class(self):
        if self.request.user.has_perm('boadicea_auth.can_risk'):
            return BwsExtendedInputSerializer
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
             defaultValue: 0.9
           - name: brca2_mut_sensitivity
             description: BRCA2 mutation sensitivity
             required: false
             type: float
             paramType: form
             defaultValue: 0.9
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
                "timestamp": datetime.datetime.now(),
                "mutation_frequency": {population: mutation_frequency},
                "mutation_sensitivity": mutation_sensitivity,
                "cancer_incidence_rates": cancer_rates,
                "pedigree_result": []
            }

            if request.user.has_perm('boadicea_auth.can_risk'):
                risk_factor_code = validated_data.get('risk_factor_code', 0)
                prs = validated_data.get('prs', None)
                if prs is not None:
                    output['prs'] = {'alpha': prs.get('alpha'), 'beta': prs.get('beta')}
                    prs = Prs(prs.get('alpha'), prs.get('beta'))
                factors = RiskFactors.decode(risk_factor_code)
                keys = list(RiskFactors.categories.keys())
                output['risk_factors'] = {keys[idx]: val for idx, val in enumerate(factors)}
            else:
                if validated_data.get('risk_factor_code', 0) > 0:
                    logger.warning('risk factor code parameter provided without the correct permissions')
                if validated_data.get('prs', 0) != 0:
                    logger.warning('polygenic risk score parameter provided without the correct permissions')
                risk_factor_code = 0
                prs = None

            try:
                warnings = PedigreeFile.validate(pf.pedigrees)
                if len(warnings) > 0:
                    output['warnings'] = warnings
            except ValidationError as e:
                logger.error(e)
                return JsonResponse(e.detail, content_type="application/json",
                                    status=status.HTTP_400_BAD_REQUEST)

            cwd = tempfile.mkdtemp(prefix=str(request.user)+"_", dir="/tmp")
            try:
                for pedi in pf.pedigrees:
                    this_pedigree = {}
                    this_pedigree["family_id"] = pedi.famid
                    this_pedigree["proband_id"] = pedi.get_target().pid

                    calcs = Predictions(pedi, mutation_frequency=mutation_frequency,
                                        mutation_sensitivity=mutation_sensitivity, cancer_rates=cancer_rates,
                                        risk_factor_code=risk_factor_code, prs=prs,
                                        cwd=cwd, request=request)
                    # Add input parameters and calculated results as attributes to 'this_pedigree'
                    self.add_attr("version", output, calcs, output)
                    self.add_attr("mutation_probabilties", this_pedigree, calcs, output)
                    self.add_attr("cancer_risks", this_pedigree, calcs, output)
                    self.add_attr("baseline_cancer_risks", this_pedigree, calcs, output)
                    self.add_attr("lifetime_cancer_risk", this_pedigree, calcs, output)
                    self.add_attr("baseline_lifetime_cancer_risk", this_pedigree, calcs, output)
                    self.add_attr("ten_yr_cancer_risk", this_pedigree, calcs, output)
                    self.add_attr("baseline_ten_yr_cancer_risk", this_pedigree, calcs, output)

                    output["pedigree_result"].append(this_pedigree)
            finally:
                shutil.rmtree(cwd)
            output_serialiser = BwsOutputSerializer(output)
            return Response(output_serialiser.data, template_name='result_tab.html')

        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)

    def add_attr(self, attr_name, this_pedigree, calcs, output):
        ''' Utility to add attribute to calculation result. '''
        try:
            this_pedigree[attr_name] = getattr(calcs, attr_name)
        except AttributeError as e:
            if 'warnings' in output:
                output['warnings'].append(attr_name+' not provided')
            else:
                output['warnings'] = [attr_name+' not provided']
            logger.debug(attr_name+' not provided :: '+str(e))
