"""
Calculate PRS from a VCF file.

© 2023 University of Cambridge
SPDX-FileCopyrightText: 2023 University of Cambridge
SPDX-License-Identifier: GPL-3.0-or-later
"""

import io
import logging
import os
from pathlib import Path
from statistics import NormalDist
import time
import traceback

from drf_spectacular.utils import extend_schema
from rest_framework import serializers, status, parsers
from rest_framework.authentication import BasicAuthentication, \
    SessionAuthentication, TokenAuthentication
from rest_framework.exceptions import NotAcceptable, ValidationError
from rest_framework.permissions import IsAuthenticated
from rest_framework.renderers import BrowsableAPIRenderer, JSONRenderer
from rest_framework.response import Response
from rest_framework.serializers import FileField
from rest_framework.views import APIView
import vcf
import vcf2prs
from vcf2prs.exception import Vcf2PrsError
from vcf2prs.prs import Prs

from bws.rest_api import RequiredAnyPermission
from bws.settings import BC_MODEL, OC_MODEL, PC_MODEL
from bws.throttles import BurstRateThrottle, EndUserIDRateThrottle, SustainedRateThrottle
from bws.serializers import PRSField


#from django.conf import settings
logger = logging.getLogger(__name__)


class Vcf2PrsInputSerializer(serializers.Serializer):
    ''' Vcf2Prs input. '''
    vcf_file = FileField(required=True, help_text=(
        "VCF genotype file. The file should be VCF format v4.0 or v4.1. It can contain "
        "both additional samples and variants not used in the PRS and are ignored."))
    sample_name = serializers.CharField(min_length=1, max_length=40, required=False,
                                        help_text="Name of the sample in the genotype file to be used to calculate the PRS")

    bc_choices = [v for k, v in (BC_MODEL['PRS_REFERENCE_FILES']['EUROPEAN'] +
                                 BC_MODEL['PRS_REFERENCE_FILES']['AFRICAN'] +
                                 BC_MODEL['PRS_REFERENCE_FILES']['EAST_ASIAN'] +
                                 BC_MODEL['PRS_REFERENCE_FILES']['SOUTH_ASIAN'])]

    oc_choices = [v for k, v in (OC_MODEL['PRS_REFERENCE_FILES']['EUROPEAN'] +
                                 OC_MODEL['PRS_REFERENCE_FILES']['AFRICAN'] +
                                 OC_MODEL['PRS_REFERENCE_FILES']['EAST_ASIAN'] +
                                 OC_MODEL['PRS_REFERENCE_FILES']['SOUTH_ASIAN'])]

    pc_choices = [v for k, v in (PC_MODEL['PRS_REFERENCE_FILES']['EUROPEAN'] +
                                 PC_MODEL['PRS_REFERENCE_FILES']['AFRICAN'] +
                                 PC_MODEL['PRS_REFERENCE_FILES']['EAST_ASIAN'] +
                                 PC_MODEL['PRS_REFERENCE_FILES']['SOUTH_ASIAN'])]
    
    bc_prs_reference_file = serializers.ChoiceField(choices=bc_choices,
                                                    help_text="Breast cancer PRS reference file", required=False)
    oc_prs_reference_file = serializers.ChoiceField(choices=oc_choices,
                                                    help_text="Ovarian cancer PRS reference file", required=False)
    pc_prs_reference_file = serializers.ChoiceField(choices=pc_choices,
                                                    help_text="Prostate cancer PRS reference file", required=False)


class Vcf2PrsOutputSerializer(serializers.Serializer):
    """ Vcf2Prs result. """
    breast_cancer_prs = PRSField(read_only=True)
    ovarian_cancer_prs = PRSField(read_only=True)
    prostate_cancer_prs = PRSField(read_only=True)


class Vcf2PrsView(APIView):
    any_perms = ['boadicea_auth.can_risk',
                 'boadicea_auth.commercial_api_breast',
                 'boadicea_auth.commercial_api_ovarian',
                 'boadicea_auth.commercial_api_prostate']   # for RequiredAnyPermission
    parser_classes = [parsers.MultiPartParser]
    renderer_classes = (JSONRenderer, BrowsableAPIRenderer, )
    serializer_class = Vcf2PrsInputSerializer
    authentication_classes = (SessionAuthentication, BasicAuthentication, TokenAuthentication, )
    permission_classes = (IsAuthenticated, RequiredAnyPermission)
    throttle_classes = (BurstRateThrottle, SustainedRateThrottle, EndUserIDRateThrottle)

    def post(self, request):
        """
        Calculate PRS from a VCF file.
        """
        start = time.time()
        serializer = self.serializer_class(data=request.data)
        if serializer.is_valid(raise_exception=True):
            validated_data = serializer.validated_data
            vcf_file = validated_data.get("vcf_file")

            moduledir = Path(vcf2prs.__file__).parent.parent
            bc_prs_ref_file = validated_data.get("bc_prs_reference_file", None)
            oc_prs_ref_file = validated_data.get("oc_prs_reference_file", None)
            pc_prs_ref_file = validated_data.get("pc_prs_reference_file", None)
            if bc_prs_ref_file is None and oc_prs_ref_file is None and pc_prs_ref_file is None:
                raise ValidationError('No breast, ovarian or prostate cancer PRS reference file provided')
            else:
                if bc_prs_ref_file is not None:
                    bc_prs_ref_file = os.path.join(moduledir, "PRSmodels_CanRisk", bc_prs_ref_file)
                if oc_prs_ref_file is not None:
                    oc_prs_ref_file = os.path.join(moduledir, "PRSmodels_CanRisk", oc_prs_ref_file)
                if pc_prs_ref_file is not None:
                    pc_prs_ref_file = os.path.join(moduledir, "PRSmodels_CanRisk", pc_prs_ref_file)

            sample_name = validated_data.get("sample_name", None)

            try:
                data = {}
                vcfStr = io.TextIOWrapper(vcf_file.file).read()
                if bc_prs_ref_file is not None:
                    data['breast_cancer_prs'] = Vcf2PrsView.get_prs(vcfStr, bc_prs_ref_file, sample_name)

                if oc_prs_ref_file is not None:
                    data['ovarian_cancer_prs'] = Vcf2PrsView.get_prs(vcfStr, oc_prs_ref_file, sample_name)

                if pc_prs_ref_file is not None:
                    data['prostate_cancer_prs'] = Vcf2PrsView.get_prs(vcfStr, pc_prs_ref_file, sample_name)

                prs_serializer = Vcf2PrsOutputSerializer(data)
                logger.info("PRS elapsed time=" + str(time.time() - start))
                return Response(prs_serializer.data)
            except Vcf2PrsError as ex:
                logger.debug(ex)
                data = {
                    'error': str(ex),
                    'samples': self.get_samples(vcfStr)
                }
                raise NotAcceptable(data)
        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)

    @staticmethod
    def get_prs(vcfStr, ref_file, sample_name):
        try:
            vcf_stream = io.StringIO(vcfStr)
            prs = Prs(prs_file=ref_file)
            prs.calculate_prs_from_vcf(vcf_stream, sample_name)
            if sample_name is None and len(prs.raw_PRS) > 1:
                raise Vcf2PrsError("Please select a sample name to use")
            
            prs.calculate_z_from_raw(prs.raw_PRS)
            alpha = prs.alpha
            zscore = prs.z_Score[0]
        finally:
            vcf_stream.close()
        return {'alpha': alpha, 'zscore': zscore, 'percent': Zscore2PercentView.get_percentage(zscore)}

    def get_samples(self, vcf_file):
        """ Get the samples in the VCF file. """
        try:
            fsock = io.StringIO(vcf_file)
            vcf_content = vcf.Reader(fsock)
            samples = vcf_content.samples
            fsock.close()
            return samples
        except Exception:
            logging.warn(traceback.format_exc())


class ZscoreInputSerializer(serializers.Serializer):
    ''' Zscore2Percent input. '''
    zscore = serializers.FloatField(required=True)


class ZscoreOutputSerializer(serializers.Serializer):
    """ PRS represented as a percentage of those with a lower PRS. """
    percent = serializers.FloatField(min_value=0, max_value=100, read_only=True)


class Zscore2PercentView(APIView):
    any_perms = ['boadicea_auth.can_risk',
                 'boadicea_auth.commercial_api_breast',
                 'boadicea_auth.commercial_api_ovarian',
                 'boadicea_auth.commercial_api_prostate']   # for RequiredAnyPermission
    renderer_classes = (JSONRenderer, BrowsableAPIRenderer, )
    serializer_class = ZscoreInputSerializer
    authentication_classes = (SessionAuthentication, BasicAuthentication, TokenAuthentication, )
    permission_classes = (IsAuthenticated, RequiredAnyPermission)
    throttle_classes = (BurstRateThrottle, SustainedRateThrottle, EndUserIDRateThrottle)

    @extend_schema(exclude=True)        # exclude from the swagger docs
    def post(self, request):
        """
        Returns PRS represented as a percentage of those with a lower PRS.
        """
        serializer = self.serializer_class(data=request.data)
        if serializer.is_valid(raise_exception=True):
            zscore = serializer.validated_data.get("zscore")
            return Response(ZscoreOutputSerializer({"percent": Zscore2PercentView.get_percentage(zscore)}).data)
        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)

    @classmethod
    def get_percentage(cls, load):
        """
        Use NormalDist to compute cumulative standard normal distribution,
        https://stackoverflow.com/questions/809362/how-to-calculate-cumulative-normal-distribution
        @param: standard normal PRS which is normally distributed in the general population with mean
        of 0 and standard deviation of 1
        @return: PRS represented as a percentage of those with a lower PRS
        """
        return NormalDist().cdf(load) * 100.0
