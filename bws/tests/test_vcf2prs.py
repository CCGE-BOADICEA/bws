from boadicea_auth.models import UserDetails
from django.contrib.auth.models import User, Permission
from django.core.urlresolvers import reverse
from django.test.testcases import TestCase
from django.test.utils import override_settings
from django.utils.encoding import force_text
from rest_framework import status
from rest_framework.authtoken.models import Token
from rest_framework.test import APIRequestFactory, APIClient
from vcf2prs import Vcf2Prs
import inspect
import json
import os


class Vcf2PrsWebServices(TestCase):
    ''' Test the risk factors webservice '''

    def setUp(self):
        ''' Create a user and set up the test client. '''
        self.factory = APIRequestFactory()
        self.client = APIClient(enforce_csrf_checks=True)
        self.user = User.objects.create_user('testuser', email='testuser@test.com',
                                             password='testing')
        # add user details
        UserDetails.objects.create(user=self.user, job_title=UserDetails.CGEN,
                                   country='UK')
        self.user.save()
        self.permission = Permission.objects.get(name='Can risk')
        self.user.user_permissions.add(self.permission)
        self.token = Token.objects.create(user=self.user)
        self.token.save()
        self.url = reverse('internal_ws:prs')
        moduledir = os.path.dirname(inspect.getfile(Vcf2Prs().__class__))
        self.vcf_file = os.path.join(moduledir, "sample_data.vcf")
        self.vcf_data = open(self.vcf_file, "r")
        self.prs_file_name = os.path.join(moduledir, "SNPs268wposition.chr_pos.txt")

    def test_prs(self):
        ''' Test POSTing to a vcf file to get a polygenic risk score. '''
        data = {'vcf_file': self.vcf_data, 'sample_name': '4'}
        self.client.credentials(HTTP_AUTHORIZATION='Token ' + self.token.key)
        response = self.client.post(self.url, data, format='multipart',
                                    HTTP_ACCEPT="application/json")
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertContains(response, "alpha")
        self.assertContains(response, "beta")

    def test_prs_v_direct(self):
        ''' Test POSTing to a vcf file to get a polygenic risk score and
        compare with the direct call to calculate a PRS. '''
        data = {'vcf_file': self.vcf_data, 'sample_name': '3'}
        self.client.credentials(HTTP_AUTHORIZATION='Token ' + self.token.key)
        response = self.client.post(self.url, data, format='multipart',
                                    HTTP_ACCEPT="application/json")
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        content = json.loads(force_text(response.content))
        prs = Vcf2Prs(prs_file_name=self.prs_file_name, geno_file_name=self.vcf_file, sample_name='3')
        load = prs.calculatePRS()
        self.assertEqual(load, content['beta'], 'web-service and direct calculation')

    def test_prs_err(self):
        ''' Test POSTing to the 400 returned without the vcf specified. '''
        data = {'sample_name': '4'}
        self.client.credentials(HTTP_AUTHORIZATION='Token ' + self.token.key)
        response = self.client.post(self.url, data, format='multipart',
                                    HTTP_ACCEPT="application/json")
        self.assertEqual(response.status_code, status.HTTP_400_BAD_REQUEST)
        self.assertEqual(response.data['vcf_file'][0], 'This field is required.')

    @override_settings(DATA_UPLOAD_MAX_MEMORY_SIZE=1)
    def test_prs_upload_limit(self):
        ''' Test POSTing to a vcf file that is too large. '''
        data = {'vcf_file': self.vcf_data, 'sample_name': '4'}
        self.client.credentials(HTTP_AUTHORIZATION='Token ' + self.token.key)
        response = self.client.post(self.url, data, format='multipart',
                                    HTTP_ACCEPT="application/json")
        self.assertEqual(response.status_code, status.HTTP_400_BAD_REQUEST)
