''' API for the BWS REST resources. '''
from rest_framework import serializers, status
from rest_framework.renderers import BrowsableAPIRenderer, JSONRenderer
from rest_framework.response import Response
from rest_framework.views import APIView
from rest_framework_xml.renderers import XMLRenderer
from django.conf import settings
from rest_framework.authentication import BasicAuthentication,\
    SessionAuthentication, TokenAuthentication
from rest_framework.permissions import IsAuthenticated
# from boadicea.perlfunc import perlreq, perl5lib, perlfunc
from boadicea.ped import PedigreeFile
from boadicea import ped
import os
import subprocess


# http://www.boriel.com/en/2007/01/21/calling-perl-from-python/
# @perlfunc
# @perlreq('Vl.pm')
# @perl5lib("/home/MINTS/tjc29/boadicea_classic/git/BOADICEA/perl/modules/vl/")
# @perl5lib("/home/MINTS/tjc29/boadicea_classic/git/BOADICEA/perl/modules/cs/")
# def vlValidateUploadedPedigreeFile(stUploadedPedigreeFileName,
#                                    rValidationRequest,
#                                    stCurrentOperation,
#                                    rMinPedigreeSize,
#                                    rMaxPedigreeSize,
#                                    ref_stTitleLeft,
#                                    ref_stTitleRight,
#                                    ref_stErrorMode,
#                                    ref_stNextPageName):
#     pass


class BwsInputSerializer(serializers.Serializer):
    ''' Boadicea result. '''
    pedigree_data = serializers.CharField()
    mut_freq = serializers.CharField()
    cancer_rates = serializers.CharField()

    def validate(self, attrs):
        mut_freq = attrs.get('mut_freq')
        cancer_rates = attrs.get('cancer_rates')

        if mut_freq not in settings.MUTATION_FREQ:
            raise serializers.ValidationError('value of mut_freq (' + mut_freq +
                                              ') is not one of the available options.')
        elif cancer_rates not in settings.CANCER_RATES:
            raise serializers.ValidationError('value of cancer_rates (' + cancer_rates +
                                              ') is not one of the available options.')
        return attrs


class BwsOutputSerializer(serializers.Serializer):
    ''' Boadicea result. '''
    mut_freq = serializers.CharField()
    cancer_rates = serializers.CharField()
    cancer_risks = serializers.CharField()
    mutation_probabilties = serializers.CharField()


class BwsView(APIView):
    renderer_classes = (XMLRenderer, JSONRenderer, BrowsableAPIRenderer, )
    serializer_class = BwsInputSerializer
    authentication_classes = (BasicAuthentication, TokenAuthentication, )
    permission_classes = (IsAuthenticated,)

    def post(self, request):
        """
        BOADICEA Web-Service (BWS)
        ---
        parameters_strategy: merge
        response_serializer: BwsOutputSerializer
        parameters:
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
           - name: brca1_mut_search_sensitivity
             description: brca1 mutation search sensitivity (only available with mut_freq=Custom)
             required: false
             type: float
             paramType: form
           - name: brca2_mut_search_sensitivity
             description: brca2 mutation search sensitivity (only available with mut_freq=Custom)
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
           - name: output_format
             description: output format
             required: true
             type: string
             paramType: form
             defaultValue: 'Percent'
             enum: ['Percent', 'Decimal']

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
            file_obj = request.FILES.get('pedigree_data')
            if file_obj is not None:
                pedigree_data = ''
                for chunk in file_obj.chunks():
                    pedigree_data += chunk.decode("utf-8")

            pf = PedigreeFile(pedigree_data)
            population = request.POST.get('mut_freq')
            cancer_rates = request.POST.get('cancer_rates')

            # mutation probablility calculation
            ped_file = pf.write_pedigree_file(file_type=ped.MUTATION_PROBS, filepath="/tmp/test_prob.ped")
            bat_file = pf.write_batch_file(ped.MUTATION_PROBS, ped_file,
                                           population=population, filepath="/tmp/test_prob.bat")
            probs = self._run(ped.MUTATION_PROBS, bat_file, cancer_rates=cancer_rates)

            # cancer risk calculation
            ped_file = pf.write_pedigree_file(file_type=ped.CANCER_RISKS, filepath="/tmp/test_risk.ped")
            bat_file = pf.write_batch_file(ped.CANCER_RISKS, ped_file,
                                           population=population, filepath="/tmp/test_risk.bat")
            risks = self._run(ped.CANCER_RISKS, bat_file, cancer_rates=cancer_rates)

#             vlValidateUploadedPedigreeFile(file_obj, 1, 'submit', 1,
#                                            275, "", "", "errorMode", "")

            output_serialiser = BwsOutputSerializer({
                    "mut_freq": request.data['mut_freq'],
                    "cancer_rates": request.data['cancer_rates'],
                    "cancer_risks": risks,
                    "mutation_probabilties": probs})

            return Response(output_serialiser.data)
        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)

    def _run(self, process_type, bat_file, cancer_rates="UK", cwd="/tmp"):
        """
        Run a process.
        """
        from subprocess import Popen, PIPE
        prog = ""
        out = ""
        if process_type == ped.MUTATION_PROBS:
            prog = os.path.join(settings.FORTRAN_HOME, "./boadicea_probs_v10.exe")
            out = "can_probs"
        else:
            prog = os.path.join(settings.FORTRAN_HOME, "./boadicea_risks_v10.exe")
            out = "can_risks"

        process = Popen([prog,
                         bat_file,  # "Sample_Pedigrees/risks_single_person.bat",
                         os.path.join(settings.FORTRAN_HOME, "Data/locus.loc"),
                         out+".stdout",
                         out+".out",
                         os.path.join(settings.FORTRAN_HOME, "Data/incidence_rates_" + cancer_rates + ".nml")],
                        cwd=cwd,
                        stdout=PIPE)

        (output, err) = process.communicate()
        try:
            exit_code = process.wait(timeout=60*4)  # timeout in seconds
            print("EXIT CODE ("+out.replace('can_', '')+"): "+str(exit_code))

            if exit_code == 0:
                with open(os.path.join(cwd, out+".out"), 'r') as myfile:
                    data = myfile.read()
                return data
            else:
                print(err)
                print(output)
        except subprocess.TimeoutExpired:
            process.terminate()
            print("we got a timeout. exiting")
