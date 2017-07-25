from django.test.testcases import TestCase
from rest_framework.test import APIRequestFactory, APIClient
from rest_framework.authtoken.models import Token
from django.contrib.auth.models import User, Permission
from boadicea_auth.models import UserDetails
from django.core.urlresolvers import reverse
import inspect
from vcf2prs import Vcf2Prs
import os
from rest_framework import status


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

    def test_prs(self):
        ''' Test POSTing to a vcf file to get a polygenic risk score. '''
        moduledir = os.path.dirname(inspect.getfile(Vcf2Prs().__class__))
        vcf_data = open(os.path.join(moduledir, "sample_data.vcf"), "r")
        data = {'vcf_file': vcf_data, 'sample_name': '4'}
        self.client.credentials(HTTP_AUTHORIZATION='Token ' + self.token.key)
        response = self.client.post(self.url, data, format='multipart',
                                    HTTP_ACCEPT="application/json")
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertContains(response, "prs")

    def test_prs_err(self):
        ''' Test POSTing to the 400 returned without the vcf specified. '''
        data = {'sample_name': '4'}
        self.client.credentials(HTTP_AUTHORIZATION='Token ' + self.token.key)
        response = self.client.post(self.url, data, format='multipart',
                                    HTTP_ACCEPT="application/json")
        self.assertEqual(response.status_code, status.HTTP_400_BAD_REQUEST)
        self.assertEqual(response.data['vcf_file'][0], 'This field is required.')
