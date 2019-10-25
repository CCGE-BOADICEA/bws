import os

from rest_framework import permissions, serializers, status
from rest_framework.authentication import BasicAuthentication, \
    SessionAuthentication, TokenAuthentication
from rest_framework.compat import coreapi, coreschema
from rest_framework.exceptions import NotAcceptable, ValidationError
from rest_framework.permissions import IsAuthenticated
from rest_framework.renderers import BrowsableAPIRenderer, JSONRenderer
from rest_framework.response import Response
from rest_framework.schemas import ManualSchema
from rest_framework.views import APIView
from vcf2prs import Vcf2Prs, Vcf2PrsError

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

logger = logging.getLogger(__name__)


class CanRiskPermission(permissions.BasePermission):
    message = 'Cancer risk factor permission not granted'

    def has_permission(self, request, view):
        return request.user.has_perm('boadicea_auth.can_risk')


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
    permission_classes = (IsAuthenticated, CanRiskPermission)
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

            moduledir = os.path.dirname(os.path.abspath(vcf2prs.__file__))
            bc_prs_ref_file = validated_data.get("bc_prs_reference_file", None)
            oc_prs_ref_file = validated_data.get("oc_prs_reference_file", None)
            if bc_prs_ref_file is None and oc_prs_ref_file is None:
                raise ValidationError('No breast or ovarian cancer PRS reference file provided')
            elif bc_prs_ref_file is not None:
                bc_prs_ref_file = os.path.join(moduledir, "PRS_files", bc_prs_ref_file)
            elif oc_prs_ref_file is not None:
                oc_prs_ref_file = os.path.join(moduledir, "PRS_files", oc_prs_ref_file)

            sample_name = validated_data.get("sample_name", None)

            try:
                if bc_prs_ref_file is not None:
                    breast_prs = Vcf2Prs(prs_file_name=bc_prs_ref_file, geno_content=vcf_file, sample_name=sample_name)
                    _raw, bc_alpha, bc_zscore = breast_prs.calculatePRS()
                else:
                    bc_alpha = 0
                    bc_zscore = 0

                if oc_prs_ref_file is not None:
                    ovarian_prs = Vcf2Prs(prs_file_name=oc_prs_ref_file, geno_content=vcf_file, sample_name=sample_name)
                    _raw, oc_alpha, oc_zscore = ovarian_prs.calculatePRS()
                else:
                    oc_alpha = 0
                    oc_zscore = 0
                data = {
                    'breast_cancer_prs': {'alpha': bc_alpha, 'zscore': bc_zscore,
                                          'percent': self.get_percentage(bc_zscore)},
                    'ovarian_cancer_prs': {'alpha': oc_alpha, 'zscore': oc_zscore,
                                           'percent': self.get_percentage(oc_zscore)}
                }
                prs_serializer = Vcf2PrsOutputSerializer(data)
                logger.info("PRS elapsed time=" + str(time.time() - start))
                return Response(prs_serializer.data)
            except Vcf2PrsError as ex:
                data = {
                    'error': ex.args[0]['Vcf2PrsError'],
                    'samples': self.get_samples(vcf_file)
                }
                raise NotAcceptable(data)
        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)

    def get_percentage(self, load):
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
