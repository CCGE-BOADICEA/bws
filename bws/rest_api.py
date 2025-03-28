'''
API for the BWS/OWS/PWS REST resources.

© 2023 University of Cambridge
SPDX-FileCopyrightText: 2023 University of Cambridge
SPDX-License-Identifier: GPL-3.0-or-later
'''
from copy import deepcopy
import datetime
import logging
import shutil
import tempfile

from django.conf import settings
from django.core.exceptions import PermissionDenied
from django.http.response import JsonResponse
from django.utils.translation import gettext_lazy as _
from drf_spectacular.utils import extend_schema
from rest_framework import status, permissions, parsers
from rest_framework.authentication import BasicAuthentication, TokenAuthentication, SessionAuthentication
from rest_framework.exceptions import ValidationError
from rest_framework.permissions import IsAuthenticated
from rest_framework.renderers import JSONRenderer, TemplateHTMLRenderer  # , BrowsableAPIRenderer
from rest_framework.response import Response
from rest_framework.views import APIView

from bws.calc.calcs import Predictions
from bws.calc.model import ModelParams
from bws.exceptions import ModelError, PedigreeError
from bws.pedigree_file import PedigreeFile, CanRiskPedigree, Prs
from bws.risk_factors.bc import BCRiskFactors
from bws.risk_factors.oc import OCRiskFactors
from bws.risk_factors.pc import PCRiskFactors
from bws.serializers import BwsInputSerializer, OutputSerializer, OwsInputSerializer, CombinedInputSerializer, \
    CombinedOutputSerializer, PwsInputSerializer
from bws.throttles import BurstRateThrottle, EndUserIDRateThrottle, SustainedRateThrottle
from bws.person import Female


logger = logging.getLogger(__name__)


class RequiredAnyPermission(permissions.BasePermission):
    """ Check that one of the permissions is met in the class variable any_perms """
    def has_permission(self, request, view):
        def test_func(user):
            for perm in view.any_perms:
                if user.has_perm(perm):
                    return True
            raise PermissionDenied()
        return test_func(request.user)


class ModelWebServiceMixin(APIView):

    parser_classes = parsers.MultiPartParser, parsers.JSONParser, parsers.FormParser
    renderer_classes = (JSONRenderer, )
    authentication_classes = (SessionAuthentication, BasicAuthentication, TokenAuthentication, )
    permission_classes = (IsAuthenticated, RequiredAnyPermission)
    throttle_classes = [BurstRateThrottle, SustainedRateThrottle, EndUserIDRateThrottle]

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

                    mname = model_settings['NAME']
                    if isinstance(pedi.get_target(), Female) and mname == "PC":
                        continue

                    if isinstance(pedi, CanRiskPedigree):
                        # for canrisk format files check if risk factors and/or prs set in the header
                        risk_factor_code = pedi.get_rfcode(mname)

                        if prs is None or len(pf.pedigrees) > 1:
                            prs = pedi.get_prs(mname)

                    this_hgt = (pedi.hgt if hasattr(pedi, 'hgt') else -1)
                    this_mdensity = (pedi.mdensity if hasattr(pedi, 'mdensity') and pedi.mdensity is not None else None)
                    if hasattr(pedi, 'ethnicity') and pedi.ethnicity is not None:
                        if params.cancer_rates != "UK":
                            raise PedigreeError(params.cancer_rates+" cancer rates with a UK ethnicity parameter ("+pedi.ons_ethnicity.get_string()+") is not valid.")
                        this_params.ethnicity = pedi.ethnicity

                    calcs = Predictions(pedi, model_params=this_params, risk_factor_code=risk_factor_code,
                                        hgt=this_hgt, mdensity=this_mdensity, prs=prs,
                                        cwd=cwd, request=request, model_settings=model_settings)
                    target = pedi.get_target()
                    # Add input parameters and calculated results as attributes to 'this_pedigree'
                    this_pedigree = {}
                    this_pedigree["family_id"] = pedi.famid
                    this_pedigree["proband_id"] = target.pid
                    this_pedigree["risk_factors"] = self.get_risk_factors(model_settings, risk_factor_code)
                    if hasattr(pedi, 'ethnicity') and pedi.ethnicity is not None:
                        this_pedigree["ethnicity"] = this_params.ethnicity.get_group()
                        this_pedigree["ons_ethnicity"] = pedi.ons_ethnicity.get_string()

                    if mname == "BC":
                        this_pedigree["risk_factors"][_('Mammographic Density')] = \
                                            this_mdensity.get_display_str() if this_mdensity is not None else "-"
                    if mname != "PC":
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
                    if int(target.age) < 50 and mname == "BC":
                        self.add_attr("ten_yr_nhs_protocol", this_pedigree, calcs, output)

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
        mname = model_settings['NAME']
        if mname == 'BC':
            rf_cls = BCRiskFactors
        elif mname == 'OC':
            rf_cls = OCRiskFactors
        elif mname == 'PC':
            rf_cls = PCRiskFactors
        else:
            raise ModelError("MODEL NOT RECOGNISED")

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


class BwsView(ModelWebServiceMixin):
    """
    Calculates the risks of breast cancer using family history, genetic and other risk factors.
    It also calculates mutation carrier probabilities in breast cancer susceptibility genes.
    """
    any_perms = ['boadicea_auth.can_risk', 'boadicea_auth.commercial_api_breast']   # for RequiredAnyPermission
    serializer_class = BwsInputSerializer

    # @profile("profile_bws.profile")
    @extend_schema(
        request=BwsInputSerializer,
        responses=OutputSerializer,
    )
    def post(self, request):
        return self.post_to_model(request, settings.BC_MODEL)


class OwsView(ModelWebServiceMixin):
    """ Ovarian Cancer Risk Model Web-Service """
    any_perms = ['boadicea_auth.can_risk', 'boadicea_auth.commercial_api_ovarian']      # for RequiredAnyPermission
    serializer_class = OwsInputSerializer

    # @profile("profile_bws.profile")
    @extend_schema(
        request=OwsInputSerializer,
        responses=OutputSerializer,
    )
    def post(self, request):
        """
        Calculates the risks of ovarian cancer using family history, genetic and other risk factors. 
        It also calculates mutation carrier probabilities in ovarian cancer susceptibility genes.
        """
        return self.post_to_model(request, settings.OC_MODEL)


class PwsView(ModelWebServiceMixin):
    """ Prostate Cancer Risk Model Web-Service """
    any_perms = ['boadicea_auth.can_risk', 'boadicea_auth.commercial_api_prostate']     # for RequiredAnyPermission
    serializer_class = PwsInputSerializer

    # @profile("profile_bws.profile")
    @extend_schema(
        request=PwsInputSerializer,
        responses=OutputSerializer,
    )
    def post(self, request):
        """
        Calculates the risks of prostate cancer using family history, genetic and other risk factors.
        """
        return self.post_to_model(request, settings.PC_MODEL)


class CombineModelResultsView(APIView):
    """
    Combine results from breast and ovarian models to produce HTML results tab.
    """
    renderer_classes = (TemplateHTMLRenderer, )
    serializer_class = CombinedInputSerializer
    authentication_classes = (SessionAuthentication, BasicAuthentication, TokenAuthentication, )
    permission_classes = (IsAuthenticated,)
    throttle_classes = (BurstRateThrottle, SustainedRateThrottle, EndUserIDRateThrottle)

    @extend_schema(exclude=True)    # exclude from the swagger docs
    def post(self, request):
        """
        Web-service to combine results from the BOADICEA and Ovarian web-services to produce
        HTML results.
        """
        serializer = self.serializer_class(data=request.data)
        if serializer.is_valid(raise_exception=True):
            validated_data = serializer.validated_data
            output_serialiser = CombinedOutputSerializer(validated_data)
            return Response(output_serialiser.data, template_name='results/tabs.html')
