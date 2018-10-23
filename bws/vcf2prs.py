import inspect
import os

from rest_framework import serializers, status
from rest_framework.authentication import BasicAuthentication, \
    SessionAuthentication, TokenAuthentication
from rest_framework.exceptions import NotAcceptable
from rest_framework.permissions import IsAuthenticated
from rest_framework.renderers import BrowsableAPIRenderer, JSONRenderer
from rest_framework.response import Response
from rest_framework.views import APIView
from rest_framework_xml.renderers import XMLRenderer
from vcf2prs import Vcf2Prs, Vcf2PrsError

from bws.serializers import FileField
from bws.risk_factors_ws import CanRiskPermission
from bws.throttles import BurstRateThrottle, EndUserIDRateThrottle, SustainedRateThrottle
import time
import logging

logger = logging.getLogger(__name__)


class Vcf2PrsInputSerializer(serializers.Serializer):
    ''' Vcf2Prs input. '''
    sample_name = serializers.CharField(min_length=1, max_length=40, required=False)
    vcf_file = FileField(required=True)
    tabix_file = FileField(required=False)


class PrsSerializer(serializers.Serializer):
    alpha = serializers.FloatField()
    beta = serializers.FloatField()


class Vcf2PrsOutputSerializer(serializers.Serializer):
    """ Vcf2Prs result. """
    breast_cancer_prs = PrsSerializer(read_only=True)
    ovarian_cancer_prs = PrsSerializer(read_only=True)


class Vcf2PrsView(APIView):
    renderer_classes = (XMLRenderer, JSONRenderer, BrowsableAPIRenderer, )
    serializer_class = Vcf2PrsInputSerializer
    authentication_classes = (SessionAuthentication, BasicAuthentication, TokenAuthentication, )
    permission_classes = (IsAuthenticated, CanRiskPermission)
    throttle_classes = (BurstRateThrottle, SustainedRateThrottle, EndUserIDRateThrottle)

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

            sample_name = validated_data.get("sample_name", None)
            moduledir = os.path.dirname(inspect.getfile(Vcf2Prs().__class__))
            prs_file_name = os.path.join(moduledir, "sample_data/BCAC_313_PRS_info_file_2018-08-06.prs")

            try:
                breast_prs = Vcf2Prs(prs_file_name=prs_file_name, geno_content=vcf_file, sample_name=sample_name)
                _raw, bc_alpha, bc_beta = breast_prs.calculatePRS()
                # NOTE:: prs_file_name to be confirmed
                ovarian_prs = Vcf2Prs(prs_file_name=prs_file_name, geno_content=vcf_file, sample_name=sample_name)
                _raw, oc_alpha, oc_beta = ovarian_prs.calculatePRS()
                data = {
                    'breast_cancer_prs': {'alpha': bc_alpha, 'beta': bc_beta},
                    'ovarian_cancer_prs': {'alpha': oc_alpha, 'beta': oc_beta}
                }
                prs_serializer = Vcf2PrsOutputSerializer(data)
                logger.info("PRS elapsed time=" + str(time.time() - start))
                return Response(prs_serializer.data)
            except Vcf2PrsError as ex:
                raise NotAcceptable(ex.args[0]['Vcf2PrsError'])
        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)
