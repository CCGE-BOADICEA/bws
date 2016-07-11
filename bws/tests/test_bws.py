import os

from django.contrib.auth.models import User
from django.core.urlresolvers import reverse
from django.test import TestCase
from rest_framework import status
from rest_framework.authentication import TokenAuthentication
from rest_framework.authtoken.models import Token
from rest_framework.test import APIRequestFactory, APIClient

from boadicea import settings


class BwsTests(TestCase):
    TEST_DATA_DIR = os.path.join(settings.PROJECT_DIR, 'tests', 'data')

    def setUp(self):
        ''' Create a user and set up the test client. '''
        self.factory = APIRequestFactory()
        self.client = APIClient(enforce_csrf_checks=True)
        self.user = User.objects.create_user('testuser', email='testuser@test.com',
                                             password='testing')
        self.user.save()
        self.token = Token.objects.create(user=self.user)
        self.token.save()

    def test_token(self):
        ''' Test token authentication. '''
        (user, _token) = TokenAuthentication().authenticate_credentials(self.token.key)
        self.assertTrue(user.is_authenticated())

    def test_token_auth_bws(self):
        ''' Test POSTing to the BWS using token authentication. '''
        url = reverse('bws')
        pedigree_data = open(os.path.join(BwsTests.TEST_DATA_DIR, "pedigree_data.txt"), "r")
        data = {'mut_freq': 'UK', 'cancer_rates': 'UK',
                'pedigree_data': pedigree_data}
        self.client.credentials(HTTP_AUTHORIZATION='Token ' + self.token.key)
        response = self.client.post(url, data, format='multipart')
        self.assertEqual(response.status_code, status.HTTP_200_OK)

    def test_force_auth_bws(self):
        ''' Test POSTing to the BWS bypassing authentication. '''
        url = reverse('bws')
        pedigree_data = open(os.path.join(BwsTests.TEST_DATA_DIR, "pedigree_data.txt"), "r")
        data = {'mut_freq': 'UK', 'cancer_rates': 'UK',
                'pedigree_data': pedigree_data}
        self.client.force_authenticate(user=self.user)
        response = self.client.post(url, data, format='multipart')
        self.assertEqual(response.status_code, status.HTTP_200_OK)
