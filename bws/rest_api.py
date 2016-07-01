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
    msg = serializers.CharField()
    pedigree_data = serializers.CharField()
    mut_freq = serializers.CharField()
    cancer_rates = serializers.CharField()


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

            output_serialiser = BwsOutputSerializer({
                    "msg": "Found data!",
                    "pedigree_data": pedigree_data,
                    "mut_freq": request.data['mut_freq'],
                    "cancer_rates": request.data['cancer_rates']})

            self._run()
            return Response(output_serialiser.data)
        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)

    def _run(self):
        from subprocess import Popen, PIPE

        process = Popen(["./boadicea_risks_v10.exe", "Sample_Pedigrees/risks_single_person.bat",
                         "Data/locus.loc", "can_risks.stdout",
                         "can_risks.out", "Data/incidence_rates_UK.nml"],
                        cwd="/home/MINTS/tjc29/boadicea_classic/git/BOADICEA/fortran10_standalone/BOADICEA_V10_ver_2_10",
                        stdout=PIPE)

        (output, err) = process.communicate()
        exit_code = process.wait()
        print(output)
        print(err)
        print(exit_code)
