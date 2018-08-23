''' API for the BWS REST resources. '''
import datetime
import logging
import shutil
import tempfile

from django.conf import settings
from rest_framework import status
from rest_framework.authentication import BasicAuthentication, TokenAuthentication, SessionAuthentication

from rest_framework.permissions import IsAuthenticated
from rest_framework.renderers import JSONRenderer, TemplateHTMLRenderer   # , BrowsableAPIRenderer
from rest_framework.response import Response
from rest_framework.views import APIView
from rest_framework_xml.renderers import XMLRenderer
from bws.pedigree import PedigreeFile, CanRiskPedigree
from bws.calcs import Predictions
from rest_framework.exceptions import NotAcceptable, ValidationError
from bws.throttles import BurstRateThrottle, EndUserIDRateThrottle, SustainedRateThrottle
from django.http.response import JsonResponse
from bws.serializers import BwsExtendedInputSerializer, BwsInputSerializer, OutputSerializer,\
    OwsExtendedInputSerializer, OwsInputSerializer
from bws.risk_factors.bc import BCRiskFactors


logger = logging.getLogger(__name__)


class Prs(object):

    def __init__(self, alpha, beta):
        ''' Polygenic risk alpha and beta values calculated from VCF file. '''
        self.alpha = alpha
        self.beta = beta


class CsrfExemptSessionAuthentication(SessionAuthentication):

    def enforce_csrf(self, request):
        return  # To not perform the csrf check previously happening


class ModelWebServiceMixin():

    def post_to_model(self, request, model_settings):
        serializer = self.serializer_class(data=request.data)
        if serializer.is_valid(raise_exception=True):
            validated_data = serializer.validated_data
            pedigree_data = validated_data.get('pedigree_data')
            pf = PedigreeFile(pedigree_data)
            population = validated_data.get('mut_freq', 'UK')
            cancer_rates = model_settings['CANCER_RATES'].get(validated_data.get('cancer_rates'))

            if population != 'Custom':
                mutation_frequency = model_settings['MUTATION_FREQUENCIES'][population]
            else:
                mutation_frequency = {}
                for gene in model_settings['GENES']:
                    try:
                        mutation_frequency[gene] = float(validated_data.get(gene.lower() + '_mut_frequency'))
                    except TypeError:
                        raise NotAcceptable("Invalid mutation frequency for " + gene + ".")

            gts = model_settings['GENETIC_TEST_SENSITIVITY']
            mutation_sensitivity = {
                k: float(validated_data.get(k.lower() + "_mut_sensitivity", gts[k]))
                for k in gts.keys()
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
                factors = BCRiskFactors.decode(risk_factor_code)
                keys = list(BCRiskFactors.categories.keys())
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

            cwd = tempfile.mkdtemp(prefix=str(request.user)+"_", dir=settings.CWD_DIR)
            try:
                for pedi in pf.pedigrees:

                    # if canrisk format file check if risk factors set in the header
                    if isinstance(pedi, CanRiskPedigree) and 'risk_factor_code' not in request.data.keys():
                        if model_settings['NAME'] == 'BC' and hasattr(pedi, 'bc_risk_factor_code'):
                            print(pedi.bc_risk_factor_code)
                            risk_factor_code = pedi.bc_risk_factor_code
                        elif model_settings['NAME'] == 'OC' and hasattr(pedi, 'oc_risk_factor_code'):
                            print(pedi.oc_risk_factor_code)
                            risk_factor_code = pedi.oc_risk_factor_code

                    this_pedigree = {}
                    this_pedigree["family_id"] = pedi.famid
                    this_pedigree["proband_id"] = pedi.get_target().pid

                    calcs = Predictions(pedi, mutation_frequency=mutation_frequency,
                                        mutation_sensitivity=mutation_sensitivity, cancer_rates=cancer_rates,
                                        risk_factor_code=risk_factor_code, prs=prs,
                                        cwd=cwd, request=request, model_settings=model_settings)
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
            except ValidationError as e:
                logger.error(e)
                return JsonResponse(e.detail, content_type="application/json",
                                    status=status.HTTP_400_BAD_REQUEST, safe=False)
            finally:
                if(not request.user.has_perm('boadicea_auth.evaluation') and
                   not request.user.has_perm('boadicea_auth.evaluation1b')):
                    shutil.rmtree(cwd)
            output_serialiser = OutputSerializer(output)
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


class BwsView(APIView, ModelWebServiceMixin):
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
        response_serializer: OutputSerializer
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
             enum: ['UK', 'UK-version-1', 'Australia', 'Canada', 'USA-white', 'Denmark', 'Finland',
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
        return self.post_to_model(request, settings.BC_MODEL)


class OwsView(APIView, ModelWebServiceMixin):
    """
    Ovarian Model Web-Service
    """
    renderer_classes = (XMLRenderer, JSONRenderer, TemplateHTMLRenderer, )
    serializer_class = OwsExtendedInputSerializer
    authentication_classes = (CsrfExemptSessionAuthentication, BasicAuthentication, TokenAuthentication, )
    permission_classes = (IsAuthenticated,)
    throttle_classes = (BurstRateThrottle, SustainedRateThrottle, EndUserIDRateThrottle)

    def get_serializer_class(self):
        if self.request.user.has_perm('boadicea_auth.can_risk'):
            return OwsExtendedInputSerializer
        return OwsInputSerializer

    # @profile("profile_bws.profile")
    def post(self, request):
        """
        Ovarian Web-Service (BWS) used to calculate the risks of ovarian cancer and the probability
        that an individual is a carrier of cancer-associated mutations in genes (BRCA1, BRCA2, RAD51D, RAD51C, BRIP1).
        As well as the individuals pedigree, the prediction model takes as input mutation frequency and sensitivity
        for each the genes and the population to use for cancer incidence rates.
        ---
        parameters_strategy: merge
        response_serializer: OutputSerializer
        parameters:
           - name: user_id
             description: unique end user ID, e.g. IP address
             type: string
             required: true
           - name: pedigree_data
             description: CanRisk pedigree data file
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
           - name: rad51d_mut_frequency
             description: RAD51D mutation frequency (only available with mut_freq=Custom)
             required: false
             type: float
             paramType: form
           - name: rad51c_mut_frequency
             description: RAD51C mutation frequency (only available with mut_freq=Custom)
             required: false
             type: float
             paramType: form
           - name: brip1_mut_frequency
             description: BRIP1 mutation frequency (only available with mut_freq=Custom)
             required: false
             type: float
             paramType: form
           - name: cancer_rates
             description: cancer incidence rates
             required: true
             type: string
             paramType: form
             defaultValue: 'UK'
             enum: ['UK', 'UK-version-1', 'Australia', 'Canada', 'USA-white', 'Denmark', 'Finland',
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
           - name: rad51d_mut_sensitivity
             description: RAD51D mutation sensitivity
             required: false
             type: float
             paramType: form
             defaultValue: 0.9
           - name: rad51c_mut_sensitivity
             description: RAD51C mutation sensitivity
             required: false
             type: float
             paramType: form
             defaultValue: 0.9
           - name: brip1_mut_sensitivity
             description: BRIP1 mutation sensitivity
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
        return self.post_to_model(request, settings.OC_MODEL)
