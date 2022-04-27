import os

from rest_framework import serializers, status
from rest_framework.authentication import BasicAuthentication, \
    SessionAuthentication, TokenAuthentication
from rest_framework.compat import coreapi, coreschema
from rest_framework.exceptions import NotAcceptable, ValidationError
from rest_framework.permissions import IsAuthenticated
from rest_framework.renderers import BrowsableAPIRenderer, JSONRenderer
from rest_framework.response import Response
from rest_framework.schemas import ManualSchema
from rest_framework.views import APIView

from bws.serializers import FileField
from bws.throttles import BurstRateThrottle, EndUserIDRateThrottle, SustainedRateThrottle
import time
import logging
import io
import vcf
import traceback
from math import erf, sqrt
from django.conf import settings
import vcf2prs
from vcf2prs.prs import Prs
from vcf2prs.exception import Vcf2PrsError
from pathlib import Path


logger = logging.getLogger(__name__)


class Vcf2PrsInputSerializer(serializers.Serializer):
    ''' Vcf2Prs input. '''
    sample_name = serializers.CharField(min_length=1, max_length=40, required=False)
    vcf_file = FileField(required=True)
    bc_prs_reference_file = serializers.CharField(min_length=1, max_length=40, required=False)
    oc_prs_reference_file = serializers.CharField(min_length=1, max_length=40, required=False)


class PrsSerializer(serializers.Serializer):
    alpha = serializers.FloatField()
    zscore = serializers.FloatField()
    percent = serializers.FloatField()


class Vcf2PrsOutputSerializer(serializers.Serializer):
    """ Vcf2Prs result. """
    breast_cancer_prs = PrsSerializer(read_only=True)
    ovarian_cancer_prs = PrsSerializer(read_only=True)


class Vcf2PrsView(APIView):
    renderer_classes = (JSONRenderer, BrowsableAPIRenderer, )
    serializer_class = Vcf2PrsInputSerializer
    authentication_classes = (SessionAuthentication, BasicAuthentication, TokenAuthentication, )
    permission_classes = (IsAuthenticated,)
    throttle_classes = (BurstRateThrottle, SustainedRateThrottle, EndUserIDRateThrottle)
    if coreapi is not None and coreschema is not None:
        schema = ManualSchema(
            fields=[
                coreapi.Field(
                    name="vcf_file",
                    required=True,
                    location='form',
                    schema=coreschema.String(
                        title="VCF",
                        description="VCF File Format",
                        format='textarea',
                    ),
                ),
                coreapi.Field(
                    name="sample_name",
                    required=True,
                    location='form',
                    schema=coreschema.String(
                        title="Sample Name",
                        description="Sample Name",
                    ),
                ),
                coreapi.Field(
                    name="bc_prs_reference_file",
                    location='form',
                    schema=coreschema.Enum(
                        list(settings.BC_MODEL['PRS_REFERENCE_FILES'].values()),
                        title="Breast cancer PRS reference file",
                        description="Breast cancer PRS reference file",
                        default=None,
                    ),
                ),
                coreapi.Field(
                    name="oc_prs_reference_file",
                    location='form',
                    schema=coreschema.Enum(
                        list(settings.OC_MODEL['PRS_REFERENCE_FILES'].values()),
                        title="Ovarian cancer PRS reference file",
                        description="Ovarian cancer PRS reference file",
                        default=None,
                    ),
                ),
            ],
            encoding="application/json",
            description="""
Variant Call Format (VCF) file to Polygenic Risk Score (PRS) web-service.
"""
        )

    def post(self, request):
        """
        Calculate PRS from a vcf file.
        ---
        response_serializer: Vcf2PrsOutputSerializer
        parameters:
           - name: sample_name
             description: name of the sample in the genotype file to be used to calculate the PRS
             type: string
             required: true
           - name: vcf_file
             description: VCF genotype file
             type: file
             required: true

        responseMessages:
           - code: 401
             message: Not authenticated

        consumes:
           - application/json
           - application/xml
        produces: ['application/json', 'application/xml']
        """
        start = time.time()
        serializer = self.serializer_class(data=request.data)
        if serializer.is_valid(raise_exception=True):
            validated_data = serializer.validated_data
            vcf_file = validated_data.get("vcf_file")
            vcf_stream = io.StringIO(vcf_file)

            moduledir = Path(vcf2prs.__file__).parent.parent
            bc_prs_ref_file = validated_data.get("bc_prs_reference_file", None)
            oc_prs_ref_file = validated_data.get("oc_prs_reference_file", None)
            if bc_prs_ref_file is None and oc_prs_ref_file is None:
                raise ValidationError('No breast or ovarian cancer PRS reference file provided')
            else:
                if bc_prs_ref_file is not None:
                    bc_prs_ref_file = os.path.join(moduledir, "PRS_files", bc_prs_ref_file)
                if oc_prs_ref_file is not None:
                    oc_prs_ref_file = os.path.join(moduledir, "PRS_files", oc_prs_ref_file)

            sample_name = validated_data.get("sample_name", None)

            try:
                if bc_prs_ref_file is not None:
                    breast_prs = Prs(prs_file=bc_prs_ref_file, geno_file=vcf_stream, sample=sample_name)
                    bc_alpha = breast_prs.alpha
                    bc_zscore = breast_prs.z_Score
                else:
                    bc_alpha = 0
                    bc_zscore = 0

                if oc_prs_ref_file is not None:
                    ovarian_prs = Prs(prs_file=oc_prs_ref_file, geno_file=vcf_stream, sample=sample_name)
                    oc_alpha = ovarian_prs.alpha
                    oc_zscore = ovarian_prs.z_Score
                else:
                    oc_alpha = 0
                    oc_zscore = 0
                data = {
                    'breast_cancer_prs': {'alpha': bc_alpha, 'zscore': bc_zscore,
                                          'percent': Zscore2PercentView.get_percentage(bc_zscore)},
                    'ovarian_cancer_prs': {'alpha': oc_alpha, 'zscore': oc_zscore,
                                           'percent': Zscore2PercentView.get_percentage(oc_zscore)}
                }
                prs_serializer = Vcf2PrsOutputSerializer(data)
                logger.info("PRS elapsed time=" + str(time.time() - start))
                return Response(prs_serializer.data)
            except Vcf2PrsError as ex:
                logger.debug(ex)
                data = {
                    'error': str(ex),
                    'samples': self.get_samples(vcf_file)
                }
                raise NotAcceptable(data)
        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)

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
    renderer_classes = (JSONRenderer, BrowsableAPIRenderer, )
    serializer_class = ZscoreInputSerializer
    authentication_classes = (SessionAuthentication, BasicAuthentication, TokenAuthentication, )
    permission_classes = (IsAuthenticated,)
    throttle_classes = (BurstRateThrottle, SustainedRateThrottle, EndUserIDRateThrottle)
    if coreapi is not None and coreschema is not None:
        schema = ManualSchema(
            fields=[
                coreapi.Field(
                    name="sample_name",
                    required=True,
                    location='form',
                    schema=coreschema.Number(
                        title="z-score",
                        description="Standard normal PRS",
                    ),
                ),
            ],
            encoding="application/json",
            description="""
Returns PRS represented as a percentage of those with a lower PRS.
"""
        )

    def post(self, request):
        """
        Returns PRS represented as a percentage of those with a lower PRS.
        ---
        response_serializer: Vcf2PrsOutputSerializer
        parameters:
           - name: z-score
             description: Standard normal PRS
             type: number
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
            zscore = serializer.validated_data.get("zscore")
            return Response(ZscoreOutputSerializer({"percent": Zscore2PercentView.get_percentage(zscore)}).data)
        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)

    @classmethod
    def get_percentage(cls, load):
        """
        Use error function to compute cumulative standard normal distribution,
        https://docs.python.org/3/library/math.html#math.erf
        (alternative to scipy.stats.norm.cdf(load, mu, sigma) * 100.0) to get
        percentage representation.
        @param: standard normal PRS which is normally distributed in the general population with mean
        of 0 and standard deviation of 1
        @return: PRS represented as a percentage of those with a lower PRS
        """
        return ((1.0 + erf(load / sqrt(2.0))) / 2.0) * 100.0
