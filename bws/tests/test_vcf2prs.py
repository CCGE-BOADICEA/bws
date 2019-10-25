from boadicea_auth.models import UserDetails
from django.contrib.auth.models import User, Permission
from django.core.urlresolvers import reverse
from django.test.testcases import TestCase
from django.test.utils import override_settings
from django.utils.encoding import force_text
from rest_framework import status
from rest_framework.authtoken.models import Token
from rest_framework.test import APIClient
from vcf2prs import Vcf2Prs
import json
import os
import vcf2prs


class Vcf2PrsWebServices(TestCase):
    ''' Test the VCF to PRS webservice '''

    @classmethod
    def setUpClass(cls):
        super(Vcf2PrsWebServices, cls).setUpClass()
        cls.client = APIClient(enforce_csrf_checks=True)
        cls.user = User.objects.create_user('testuser', email='testuser@test.com',
                                            password='testing')
        # add user details
        UserDetails.objects.create(user=cls.user, job_title=UserDetails.CGEN,
                                   country='UK')
        cls.user.save()
        cls.permission = Permission.objects.get(name='Can risk')
        cls.user.user_permissions.add(cls.permission)
        cls.token = Token.objects.create(user=cls.user)
        cls.token.save()
        cls.url = reverse('prs')
        cls.moduledir = os.path.dirname(os.path.abspath(vcf2prs.__file__))

    def setUp(self):
        ''' Create a user and set up the test client. '''
        self.vcf_file = os.path.join(Vcf2PrsWebServices.moduledir, "sample_VCF_files/sample_BCAC_313.vcf")
        self.vcf_data = open(self.vcf_file, "r")
        self.prs_reference_file = "BCAC_313_PRS.prs"
        self.prs_file_name = os.path.join(Vcf2PrsWebServices.moduledir, "PRS_files", self.prs_reference_file)

    def test_prs(self):
        ''' Test POSTing to a vcf file to get a polygenic risk score. '''
        data = {'vcf_file': self.vcf_data, 'sample_name': 'Mod', 'bc_prs_reference_file': self.prs_reference_file}
        Vcf2PrsWebServices.client.credentials(HTTP_AUTHORIZATION='Token ' + Vcf2PrsWebServices.token.key)
        response = Vcf2PrsWebServices.client.post(Vcf2PrsWebServices.url, data, format='multipart',
                                                  HTTP_ACCEPT="application/json")
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertContains(response, "alpha")
        self.assertContains(response, "zscore")

    def test_prs_v_direct(self):
        ''' Test POSTing to a vcf file to get a polygenic risk score and
        compare with the direct call to calculate a PRS. '''
        data = {'vcf_file': self.vcf_data, 'sample_name': 'High', 'bc_prs_reference_file': self.prs_reference_file}
        Vcf2PrsWebServices.client.credentials(HTTP_AUTHORIZATION='Token ' + Vcf2PrsWebServices.token.key)
        response = Vcf2PrsWebServices.client.post(Vcf2PrsWebServices.url, data, format='multipart',
                                                  HTTP_ACCEPT="application/json")
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        content = json.loads(force_text(response.content))
        prs = Vcf2Prs(prs_file_name=self.prs_file_name, geno_file_name=self.vcf_file, sample_name='High')
        _raw, _alpha, zscore = prs.calculatePRS()
        self.assertEqual(zscore, content['breast_cancer_prs']['zscore'], 'web-service and direct calculation')

    def test_prs_err(self):
        ''' Test POSTing to the 400 returned without the vcf specified. '''
        data = {'sample_name': 'SampleA'}
        Vcf2PrsWebServices.client.credentials(HTTP_AUTHORIZATION='Token ' + Vcf2PrsWebServices.token.key)
        response = Vcf2PrsWebServices.client.post(Vcf2PrsWebServices.url, data, format='multipart',
                                                  HTTP_ACCEPT="application/json")
        self.assertEqual(response.status_code, status.HTTP_400_BAD_REQUEST)
        self.assertEqual(response.data['vcf_file'][0], 'This field is required.')

    @override_settings(DATA_UPLOAD_MAX_MEMORY_SIZE=1)
    def test_prs_upload_limit(self):
        ''' Test POSTing to a vcf file that is too large. '''
        data = {'vcf_file': self.vcf_data, 'sample_name': 'SampleA'}
        Vcf2PrsWebServices.client.credentials(HTTP_AUTHORIZATION='Token ' + Vcf2PrsWebServices.token.key)
        response = Vcf2PrsWebServices.client.post(Vcf2PrsWebServices.url, data, format='multipart',
                                                  HTTP_ACCEPT="application/json")
        self.assertEqual(response.status_code, status.HTTP_400_BAD_REQUEST)
