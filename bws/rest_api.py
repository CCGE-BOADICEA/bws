''' API for the BWS/OWS REST resources. '''
from copy import deepcopy
import datetime
import logging
import re
import shutil
import tempfile

from django.conf import settings
from django.http.response import JsonResponse
from django.utils.translation import gettext_lazy as _
from rest_framework import status
from rest_framework.authentication import BasicAuthentication, TokenAuthentication, SessionAuthentication
from rest_framework.compat import coreapi, coreschema
from rest_framework.exceptions import ValidationError
from rest_framework.permissions import IsAuthenticated
from rest_framework.renderers import JSONRenderer, TemplateHTMLRenderer  # , BrowsableAPIRenderer
from rest_framework.response import Response
from rest_framework.schemas import ManualSchema
from rest_framework.views import APIView

from bws.calcs import Predictions, ModelParams, Risk, ModelOpts
from bws.pedigree import PedigreeFile, CanRiskPedigree, Prs
from bws.risk_factors.bc import BCRiskFactors
from bws.risk_factors.oc import OCRiskFactors
from bws.serializers import BwsInputSerializer, OutputSerializer, OwsInputSerializer, CombinedInputSerializer, \
    CombinedOutputSerializer, BCTenYrSerializer
from bws.throttles import BurstRateThrottle, EndUserIDRateThrottle, SustainedRateThrottle


logger = logging.getLogger(__name__)


class ModelWebServiceMixin():

    def post_to_model(self, request, model_settings):
        serializer = self.serializer_class(data=request.data)
        if serializer.is_valid(raise_exception=True):
            validated_data = serializer.validated_data
            pf = PedigreeFile(validated_data.get('pedigree_data'))
            params = ModelParams.factory(validated_data, model_settings)

            output = {
                "timestamp": datetime.datetime.now(),
                "mutation_frequency": {params.population: params.mutation_frequency},
                "mutation_sensitivity": params.mutation_sensitivity,
                "cancer_incidence_rates": params.cancer_rates,
                "pedigree_result": []
            }

            prs = validated_data.get('prs', None)
            if prs is not None:
                prs = Prs(prs.get('alpha'), prs.get('zscore'))

            try:
                warnings = PedigreeFile.validate(pf.pedigrees)
                if len(warnings) > 0:
                    output['warnings'] = warnings
            except ValidationError as e:
                logger.error(e)
                return JsonResponse(e.detail, content_type="application/json", status=status.HTTP_400_BAD_REQUEST)

            # note limit username string length used here to avoid paths too long for model code
            cwd = tempfile.mkdtemp(prefix=str(request.user)[:20]+"_", dir=settings.CWD_DIR)
            try:
                for pedi in pf.pedigrees:
                    risk_factor_code = 0
                    this_params = deepcopy(params)
                    # check if Ashkenazi Jewish status set & correct mutation frequencies
                    if pedi.is_ashkn() and not settings.REGEX_ASHKN.match(params.population):
                        msg = 'mutation frequencies set to Ashkenazi Jewish population values ' \
                              'for family ('+pedi.famid+') as a family member has Ashkenazi Jewish status.'
                        logger.debug('mutation frequencies set to Ashkenazi Jewish population values')
                        if 'warnings' in output:
                            output['warnings'].append(msg)
                        else:
                            output['warnings'] = [msg]
                        this_params.isashk = True
                        this_params.population = 'Ashkenazi'
                        this_params.mutation_frequency = model_settings['MUTATION_FREQUENCIES']['Ashkenazi']

                    if isinstance(pedi, CanRiskPedigree):
                        # for canrisk format files check if risk factors and/or prs set in the header
                        mname = model_settings['NAME']
                        risk_factor_code = pedi.get_rfcode(mname)

                        if prs is None or len(pf.pedigrees) > 1:
                            prs = pedi.get_prs(mname)

                    this_hgt = (pedi.hgt if hasattr(pedi, 'hgt') else -1)
                    calcs = Predictions(pedi, model_params=this_params,
                                        risk_factor_code=risk_factor_code, hgt=this_hgt, prs=prs,
                                        cwd=cwd, request=request, model_settings=model_settings)
                    # Add input parameters and calculated results as attributes to 'this_pedigree'
                    this_pedigree = {}
                    this_pedigree["family_id"] = pedi.famid
                    this_pedigree["proband_id"] = pedi.get_target().pid
                    this_pedigree["risk_factors"] = self.get_risk_factors(model_settings, risk_factor_code)
                    this_pedigree["risk_factors"][_('Height (cm)')] = this_hgt if this_hgt != -1 else "-"
                    if prs is not None:
                        this_pedigree["prs"] = {'alpha': prs.alpha, 'zscore': prs.zscore}
                    this_pedigree["mutation_frequency"] = {this_params.population: this_params.mutation_frequency}
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
                shutil.rmtree(cwd)
                # print(model_settings['NAME']+" :: "+cwd)
            output_serialiser = OutputSerializer(output)
            return Response(output_serialiser.data, template_name='result_tab_gp.html')

        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)

    def get_risk_factors(self, model_settings, risk_factor_code):
        ''' Get a dictionary of the decoded risk factor categories from the risk factor code. '''
        rf_cls = BCRiskFactors if model_settings['NAME'] == 'BC' else OCRiskFactors
        rfcats = rf_cls.decode(risk_factor_code)
        rfs = rf_cls.risk_factors
        return {rfs[idx].snake_name(): rfs[idx].cats[val] for idx, val in enumerate(rfcats)}

    def add_attr(self, attr_name, this_pedigree, calcs, output):
        ''' Utility to add attribute to calculation result. '''
        try:
            this_pedigree[attr_name] = getattr(calcs, attr_name)
        except AttributeError as e:
            if 'warnings' in output:
                output['warnings'].append(attr_name+' not provided')
            else:
                output['warnings'] = [attr_name+' not provided']
            logger.debug(f'{attr_name} not provided :: {e}')

    @classmethod
    def get_fields(cls, model):
        ''' Generate schema fields used to generate API docs. '''
        fields = [
            coreapi.Field(
                name="pedigree_data",
                required=True,
                location='form',
                schema=coreschema.String(
                    title="Pedigree",
                    description="CanRisk File Format",
                    format='textarea',
                ),
            ),
            coreapi.Field(
                name="user_id",
                required=True,
                location='form',
                schema=coreschema.String(
                    title="User ID",
                    description="Unique end user ID",
                ),
            ),
            coreapi.Field(
                name="cancer_rates",
                required=True,
                location='form',
                schema=coreschema.Enum(
                    list(settings.BC_MODEL['CANCER_RATES'].keys()),
                    title="Cancer rates",
                    description="Cancer incidence rates",
                    default="UK",
                ),
            ),
            coreapi.Field(
                name="mut_freq",
                required=True,
                location='form',
                schema=coreschema.Enum(
                    list(settings.BC_MODEL['MUTATION_FREQUENCIES'].keys()),
                    title="Mutation frequency",
                    description="Mutation frequency",
                    default="UK",
                ),
            ),
            coreapi.Field(
                name="prs",
                required=False,
                location='form',
                schema=coreschema.Object(
                    title="Polygenic risk score",
                    description='PRS, e.g. {"alpha":0.45,"zscore":2.652}',
                    properties={'alpha': coreschema.Number, 'zscore': coreschema.Number},
                ),
            ),
            # coreapi.Field(
            #    name="risk_factor_code",
            #    required=False,
            #    location='form',
            #    schema=coreschema.Integer(
            #        minimum=0,
            #        description="Risk factor code",
            #    ),
            # ),
        ]

        # fields += [
        #    coreapi.Field(
        #        name=g.lower() + "_mut_frequency",
        #        required=False,
        #        location='form',
        #        schema=coreschema.Number(
        #            title=g+" mutation frequency",
        #            description=g+' mutation frequency',
        #            minimum=settings.MIN_MUTATION_FREQ,
        #            maximum=settings.MAX_MUTATION_FREQ
        #        ),
        #    ) for g in model['GENES']
        # ]
        fields += [
            coreapi.Field(
                name=g.lower() + "_mut_sensitivity",
                required=False,
                location='form',
                schema=coreschema.Number(
                    title=g+" mutation sensitivity",
                    description=g+' mutation sensitivity',
                    maximum=1,
                ),
            ) for g in model['GENES']
        ]
        return fields


class BwsView(APIView, ModelWebServiceMixin):
    """
    BOADICEA Web-Service
    """
    renderer_classes = (JSONRenderer, TemplateHTMLRenderer, )
    serializer_class = BwsInputSerializer
    authentication_classes = (SessionAuthentication, BasicAuthentication, TokenAuthentication, )
    permission_classes = (IsAuthenticated,)
    throttle_classes = (BurstRateThrottle, SustainedRateThrottle, EndUserIDRateThrottle)
    model = settings.BC_MODEL
    if coreapi is not None and coreschema is not None:
        schema = ManualSchema(
            fields=ModelWebServiceMixin.get_fields(model),
            encoding="application/json",
            description="""
BOADICEA Web-Service (BWS) used to calculate the risks of breast cancer and the probability
that an individual is a carrier of cancer-associated mutations in genes (""" + ', '.join([g for g in model['GENES']]) + """).
As well as the individuals pedigree, the prediction model takes as input mutation frequency and sensitivity
for each the genes and the population to use for cancer incidence rates.
"""
        )

    # @profile("profile_bws.profile")
    def post(self, request):
        """
        BOADICEA Web-Service (BWS) used to calculate the risks of breast cancer and the probability
        that an individual is a carrier of cancer-associated mutations in genes (BRCA1, BRCA2, PALB2, CHEK2, ATM...).
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
             enum: ['UK', 'Ashkenazi', 'Iceland']
           - name: cancer_rates
             description: cancer incidence rates
             required: true
             type: string
             paramType: form
             defaultValue: 'UK'
             enum: ['UK', 'Australia', 'Canada', 'USA', 'Denmark', 'Estonia', 'Finland', 'France',
                    'Iceland', 'Netherlands', 'New-Zealand', 'Norway', 'Slovenia', 'Spain', 'Sweden']
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
    renderer_classes = (JSONRenderer, TemplateHTMLRenderer, )
    serializer_class = OwsInputSerializer
    authentication_classes = (SessionAuthentication, BasicAuthentication, TokenAuthentication, )
    permission_classes = (IsAuthenticated,)
    throttle_classes = (BurstRateThrottle, SustainedRateThrottle, EndUserIDRateThrottle)
    model = settings.OC_MODEL
    if coreapi is not None and coreschema is not None:
        schema = ManualSchema(
                fields=ModelWebServiceMixin.get_fields(model),
                encoding="application/json",
                description="""
Ovarian Web-Service (OWS) used to calculate the risks of ovarian cancer and the probability
that an individual is a carrier of cancer-associated mutations in genes (""" + ', '.join([g for g in model['GENES']]) + """).
As well as the individuals pedigree, the prediction model takes as input mutation frequency and sensitivity
for each the genes and the population to use for cancer incidence rates.
"""
            )

    # @profile("profile_bws.profile")
    def post(self, request):
        """
        Ovarian Web-Service (OWS) used to calculate the risks of ovarian cancer and the probability
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
             enum: ['UK', 'Ashkenazi', 'Iceland']
           - name: cancer_rates
             description: cancer incidence rates
             required: true
             type: string
             paramType: form
             defaultValue: 'UK'
             enum: ['UK', 'Australia', 'Canada', 'USA', 'Denmark', 'Estonia', 'Finland', 'France',
                    'Iceland', 'Netherlands', 'New-Zealand', 'Norway', 'Slovenia', 'Spain', 'Sweden']
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


class BCTenYrView(APIView, ModelWebServiceMixin):
    """
    Ten year breast cancer risks calculation Web-Service
    """
    renderer_classes = (JSONRenderer, TemplateHTMLRenderer, )
    serializer_class = BCTenYrSerializer
    authentication_classes = (SessionAuthentication, BasicAuthentication, TokenAuthentication, )
    permission_classes = (IsAuthenticated,)
    throttle_classes = (BurstRateThrottle, SustainedRateThrottle, EndUserIDRateThrottle)
    model = settings.BC_MODEL
    if coreapi is not None and coreschema is not None:
        fields = ModelWebServiceMixin.get_fields(model)
        fields.insert(0, coreapi.Field(
                name="tenyr_ages",
                required=True,
                location='form',
                schema=coreschema.Array(
                    title="tenyr_ages",
                    description='List of ages to calculate the 10-year risks',
                    min_items=1,
                    unique_items=True)
            )
        )
        schema = ManualSchema(
            fields=fields,
            encoding="application/json",
            description="""
Ten year breast cancer risks calculations as per those given for the ages 40-49
(https://canrisk.atlassian.net/wiki/x/NwDCAg). The web-service takes a list of
ages to calculate the 10-year risks for, e.g. [25, 26, 27, 28, 29] or [29].
"""
        )

    # @profile("profile_bws.profile")
    def post(self, request):
        serializer = self.serializer_class(data=request.data)
        if serializer.is_valid(raise_exception=True):
            validated_data = serializer.validated_data
            pf = PedigreeFile(validated_data.get('pedigree_data'))
            model_settings = settings.BC_MODEL
            params = ModelParams.factory(validated_data, model_settings)

            tenyr_ages = re.sub("[\[\]]", "", validated_data.get('tenyr_ages'))
            tenyr_ages = [int(item.strip()) for item in tenyr_ages.split(',')]

            output = {
                "timestamp": datetime.datetime.now(),
                "mutation_frequency": {params.population: params.mutation_frequency},
                "mutation_sensitivity": params.mutation_sensitivity,
                "cancer_incidence_rates": params.cancer_rates,
                "pedigree_result": []
            }

            prs = validated_data.get('prs', None)
            if prs is not None:
                prs = Prs(prs.get('alpha'), prs.get('zscore'))

            try:
                warnings = PedigreeFile.validate(pf.pedigrees)
                if len(warnings) > 0:
                    output['warnings'] = warnings
            except ValidationError as e:
                logger.error(e)
                return JsonResponse(e.detail, content_type="application/json", status=status.HTTP_400_BAD_REQUEST)

            # note limit username string length used here to avoid paths too long for model code
            cwd = tempfile.mkdtemp(prefix=str(request.user)[:20]+"_", dir=settings.CWD_DIR)
            try:
                for pedi in pf.pedigrees:
                    risk_factor_code = 0
                    this_params = deepcopy(params)
                    # check if Ashkenazi Jewish status set & correct mutation frequencies
                    if pedi.is_ashkn() and not settings.REGEX_ASHKN.match(params.population):
                        msg = 'mutation frequencies set to Ashkenazi Jewish population values ' \
                              'for family ('+pedi.famid+') as a family member has Ashkenazi Jewish status.'
                        logger.debug('mutation frequencies set to Ashkenazi Jewish population values')
                        if 'warnings' in output:
                            output['warnings'].append(msg)
                        else:
                            output['warnings'] = [msg]
                        this_params.isashk = True
                        this_params.population = 'Ashkenazi'
                        this_params.mutation_frequency = model_settings['MUTATION_FREQUENCIES']['Ashkenazi']

                    if isinstance(pedi, CanRiskPedigree):
                        # for canrisk format files check if risk factors and/or prs set in the header
                        mname = model_settings['NAME']
                        risk_factor_code = pedi.get_rfcode(mname)

                        if prs is None or len(pf.pedigrees) > 1:
                            prs = pedi.get_prs(mname)

                    this_hgt = (pedi.hgt if hasattr(pedi, 'hgt') else -1)
                    calcs = Predictions(pedi, model_params=this_params,
                                        risk_factor_code=risk_factor_code, hgt=this_hgt, prs=prs, run_risks=False,
                                        cwd=cwd, request=request, model_settings=model_settings)
                    calcs.niceness = Predictions._get_niceness(calcs.pedi)

                    calcs.ten_yr_cancer_risk = []
                    t = calcs.pedi.get_target()
                    if not t.cancers.is_cancer_diagnosed():
                        model_opts = ModelOpts(out="rr_10yr.txt", probs=False, rj=True, rl=False, rr=False, ry=False)
                        _rl, _rr, _ry, rj, _mp = Risk(calcs).get_risk(model_opts)
                        if rj is not None:
                            calcs.ten_yr_cancer_risk = rj

                    # Add input parameters and calculated results as attributes to 'this_pedigree'
                    this_pedigree = {}
                    this_pedigree["family_id"] = pedi.famid
                    this_pedigree["proband_id"] = pedi.get_target().pid
                    this_pedigree["risk_factors"] = self.get_risk_factors(model_settings, risk_factor_code)
                    this_pedigree["risk_factors"][_('Height (cm)')] = this_hgt if this_hgt != -1 else "-"
                    if prs is not None:
                        this_pedigree["prs"] = {'alpha': prs.alpha, 'zscore': prs.zscore}
                    this_pedigree["mutation_frequency"] = {this_params.population: this_params.mutation_frequency}
                    self.add_attr("version", output, calcs, output)
                    self.add_attr("ten_yr_cancer_risk", this_pedigree, calcs, output)

                    output["pedigree_result"].append(this_pedigree)
            except ValidationError as e:
                logger.error(e)
                return JsonResponse(e.detail, content_type="application/json",
                                    status=status.HTTP_400_BAD_REQUEST, safe=False)
            finally:
                shutil.rmtree(cwd)
                # print("BCTenYr :: "+cwd)
            output_serialiser = OutputSerializer(output)
            return Response(output_serialiser.data, template_name='result_tab_gp.html')

        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)


class CombineModelResultsView(APIView):
    """
    Combine results from breast and ovarian models to produce HTML results tab.
    """
    renderer_classes = (TemplateHTMLRenderer, )
    serializer_class = CombinedInputSerializer
    authentication_classes = (SessionAuthentication, BasicAuthentication, TokenAuthentication, )
    permission_classes = (IsAuthenticated,)
    throttle_classes = (BurstRateThrottle, SustainedRateThrottle, EndUserIDRateThrottle)

    def post(self, request):
        """
        Web-service to combine results from the BOADICEA and Ovarian web-services to produce
        HTML results.
        ---
        parameters_strategy: merge
        response_serializer: OutputSerializer
        parameters:
           - name: ows_result
             description: ovarian cancer web service result
             ptype: OutputSerializer
             required: true
           - name: bws_result
             description: breast cancer web service result
             ptype: OutputSerializer
             required: true

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
            output_serialiser = CombinedOutputSerializer(validated_data)
            return Response(output_serialiser.data, template_name='result_tab.html')
