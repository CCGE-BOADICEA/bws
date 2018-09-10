""" BOADICEA web-service testing.  """
import os

from django.conf import settings
from django.contrib.auth.models import User, Permission
from django.core.urlresolvers import reverse
from django.test import TestCase
from rest_framework import status
from rest_framework.authtoken.models import Token
from rest_framework.test import APIClient
from django.utils.encoding import force_text
import json
from boadicea_auth.models import UserDetails
from django.test.utils import override_settings
from bws.risk_factors.bc import BCRiskFactors


class BwsTests(TestCase):
    TEST_BASE_DIR = os.path.dirname(os.path.dirname(__file__))
    TEST_DATA_DIR = os.path.join(TEST_BASE_DIR, 'tests', 'data')

    @classmethod
    def setUpClass(cls):
        ''' Create a user and set up the test client. '''
        super(BwsTests, cls).setUpClass()
        cls.user = User.objects.create_user('testuser', email='testuser@test.com',
                                            password='testing')
        # add user details
        UserDetails.objects.create(user=cls.user, job_title=UserDetails.CGEN,
                                   country='UK')
        cls.user.save()
        cls.token = Token.objects.create(user=cls.user)
        cls.token.save()
        cls.url = reverse('bws')

    def setUp(self):
        ''' Set up pedigree data. '''
        self.client = APIClient(enforce_csrf_checks=True)
        self.pedigree_data = open(os.path.join(BwsTests.TEST_DATA_DIR, "pedigree_data.txt"), "r")

    def test_token_auth_bws(self):
        ''' Test POSTing to the BWS using token authentication. '''
        data = {'mut_freq': 'UK', 'cancer_rates': 'UK',
                'pedigree_data': self.pedigree_data,
                'user_id': 'test_XXX'}
        self.client.credentials(HTTP_AUTHORIZATION='Token ' + BwsTests.token.key)
        response = self.client.post(BwsTests.url, data, format='multipart',
                                    HTTP_ACCEPT="application/json")
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        content = json.loads(force_text(response.content))
        self.assertTrue("mutation_frequency" in content)
        self.assertTrue("pedigree_result" in content)
        self.assertTrue("family_id" in content["pedigree_result"][0])

    def test_multi_pedigree_bws(self):
        ''' Test POSTing multiple pedigrees to the BWS. '''
        multi_pedigree_data = open(os.path.join(BwsTests.TEST_DATA_DIR, "multi_pedigree_data.txt"), "r")
        data = {'mut_freq': 'UK', 'cancer_rates': 'UK',
                'pedigree_data': multi_pedigree_data,
                'user_id': 'test_XXX'}
        self.client.credentials(HTTP_AUTHORIZATION='Token ' + BwsTests.token.key)
        response = self.client.post(BwsTests.url, data, format='multipart',
                                    HTTP_ACCEPT="application/json")
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        content = json.loads(force_text(response.content))
        self.assertEqual(len(content['pedigree_result']), 2, "two results")
        family_ids = ["XXX0", "XXX1"]
        for res in content['pedigree_result']:
            self.assertTrue(res['family_id'] in family_ids)

    def test_canrisk_format_bws(self):
        ''' Test POSTing canrisk format pedigree to the BWS. '''
        canrisk_data = open(os.path.join(BwsTests.TEST_DATA_DIR, "canrisk_data_v1.txt"), "r")
        data = {'mut_freq': 'UK', 'cancer_rates': 'UK',
                'pedigree_data': canrisk_data,
                'user_id': 'test_XXX'}
        self.client.credentials(HTTP_AUTHORIZATION='Token ' + BwsTests.token.key)
        response = self.client.post(BwsTests.url, data, format='multipart',
                                    HTTP_ACCEPT="application/json")
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        content = json.loads(force_text(response.content))
        self.assertTrue("pedigree_result" in content)
        genes = settings.BC_MODEL['GENES']
        for g in genes:
            self.assertTrue(g in content['mutation_frequency']['UK'])

#     def test_token_auth_err(self):
#         ''' Test POSTing to the BWS using token authentication. '''
#         data = {'mut_freq': 'UK', 'cancer_rates': 'UK',
#                 'pedigree_data': self.pedigree_data}
#         self.client.credentials(HTTP_AUTHORIZATION='Token xxxxxxxxxx')
#         response = self.client.post(BwsTests.url, data, format='multipart',
#                                         HTTP_ACCEPT="application/json")
#         self.assertNotEqual(response.status_code, status.HTTP_200_OK)
#         content = json.loads(force_text(response.content))
#         self.assertEqual(content['detail'], 'Invalid token.')

    def test_force_auth_bws(self):
        ''' Test POSTing to the BWS bypassing authentication. '''
        data = {'mut_freq': 'UK', 'cancer_rates': 'UK',
                'pedigree_data': self.pedigree_data,
                'user_id': 'test_XXX'}
        self.client.force_authenticate(user=BwsTests.user)
        response = self.client.post(BwsTests.url, data, format='multipart')
        self.assertEqual(response.status_code, status.HTTP_200_OK)

    def test_custom_mutation_frequency(self):
        ''' Test POSTing custom mutation frequencies. '''
        genes = settings.BC_MODEL['GENES']
        data = {g.lower() + '_mut_frequency': 0.00085 for g in genes}
        data.update({'mut_freq': 'Custom', 'cancer_rates': 'UK',
                     'pedigree_data': self.pedigree_data, 'user_id': 'test_XXX'})

        self.client.force_authenticate(user=BwsTests.user)
        response = self.client.post(BwsTests.url, data, format='multipart',
                                    HTTP_ACCEPT="application/json")
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        content = json.loads(force_text(response.content))

        for g, mf in content['mutation_frequency']['Custom'].items():
            self.assertTrue(g in genes)
            self.assertEqual(mf, data[g.lower() + '_mut_frequency'])

    def test_custom_mutation_frequency_errs(self):
        ''' Test POSTing custom mutation frequencies with errors. '''
        genes = settings.BC_MODEL['GENES']
        data = {g.lower() + '_mut_frequency':
                (settings.MAX_MUTATION_FREQ + 0.1) if idx % 2 == 0 else (settings.MAX_MUTATION_FREQ - 0.1)
                for idx, g in enumerate(genes)}
        data.update({'mut_freq': 'Custom', 'cancer_rates': 'UK',
                     'pedigree_data': self.pedigree_data, 'user_id': 'test_XXX'})

        self.client.force_authenticate(user=BwsTests.user)
        response = self.client.post(BwsTests.url, data, format='multipart',
                                    HTTP_ACCEPT="application/json")
        self.assertEqual(response.status_code, status.HTTP_400_BAD_REQUEST)
        content = json.loads(force_text(response.content))
        self.assertEqual(len(content.keys()), len(genes))
        for k in content.keys():
            self.assertTrue(k.split("_")[0].upper() in genes)

    def test_missing_fields(self):
        ''' Test POSTing with missing required fields. '''
        data = {'mut_freq': 'UK'}
        self.client.credentials(HTTP_AUTHORIZATION='Token ' + BwsTests.token.key)
        response = self.client.post(BwsTests.url, data, format='multipart',
                                    HTTP_ACCEPT="application/json")
        self.assertEqual(response.status_code, status.HTTP_400_BAD_REQUEST)
        content = json.loads(force_text(response.content))
        self.assertEqual(content['user_id'][0], 'This field is required.')
        self.assertEqual(content['cancer_rates'][0], 'This field is required.')
        self.assertEqual(content['pedigree_data'][0], 'This field is required.')

    def test_bws_risk_factor(self):
        ''' Test affect of including the risk factors. '''
        data = {'mut_freq': 'UK', 'cancer_rates': 'UK',
                'pedigree_data': self.pedigree_data,
                'user_id': 'test_XXX', 'risk_factor_code': 7}
        self.client.credentials(HTTP_AUTHORIZATION='Token ' + BwsTests.token.key)
        # no permissions to use the risk factors and so ignored
        response = self.client.post(BwsTests.url, data, format='multipart', HTTP_ACCEPT="application/json")
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        cancer_risks1 = json.loads(force_text(response.content))['pedigree_result'][0]['cancer_risks']

        # add permissions to enable use of the risk factors
        data['pedigree_data'] = open(os.path.join(BwsTests.TEST_DATA_DIR, "pedigree_data.txt"), "r")
        BwsTests.user.user_permissions.add(Permission.objects.get(name='Can risk'))
        response = self.client.post(BwsTests.url, data, format='multipart', HTTP_ACCEPT="application/json")
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        cancer_risks2 = json.loads(force_text(response.content))['pedigree_result'][0]['cancer_risks']
        self.assertLess(cancer_risks2[0]['breast cancer risk']['decimal'],
                        cancer_risks1[0]['breast cancer risk']['decimal'])

    def test_risk_factors_inconsistent(self):
        ''' Test inconsistent risk factors, e.g. age of first birth specified with parity unobserved. '''
        data = {'mut_freq': 'UK', 'cancer_rates': 'UK',
                'pedigree_data': self.pedigree_data,
                'user_id': 'test_XXX', 'risk_factor_code': BCRiskFactors.encode([0, 0, 1, 0, 0, 0, 0, 0, 0, 0])}
        self.client.credentials(HTTP_AUTHORIZATION='Token ' + BwsTests.token.key)
        BwsTests.user.user_permissions.add(Permission.objects.get(name='Can risk'))
        response = self.client.post(BwsTests.url, data, format='multipart', HTTP_ACCEPT="application/json")
        self.assertEqual(response.status_code, status.HTTP_400_BAD_REQUEST)
        content = json.loads(force_text(response.content))
        self.assertTrue('Model Error' in content)

    def test_bws_errors(self):
        ''' Test an error is reported by the web-service for an invalid year of birth. '''
        # force an error changing to an invalid year of birth
        ped = open(os.path.join(BwsTests.TEST_DATA_DIR, "pedigree_data.txt"), "r")
        pd = ped.read().replace('1963', '1600')
        data = {'mut_freq': 'UK', 'cancer_rates': 'UK',
                'pedigree_data': pd, 'user_id': 'test_XXX'}
        self.client.force_authenticate(user=BwsTests.user)
        response = self.client.post(BwsTests.url, data, format='multipart',
                                    HTTP_ACCEPT="application/json")
        self.assertEqual(response.status_code, status.HTTP_400_BAD_REQUEST)
        content = json.loads(force_text(response.content))
        self.assertTrue('Person Error' in content)
        self.assertTrue('year of birth' in content['Person Error'])

    @override_settings(FORTRAN_TIMEOUT=0.05)
    def test_bws_timeout(self):
        ''' Test a timeout error is reported by the web-service. '''
        data = {'mut_freq': 'UK', 'cancer_rates': 'UK',
                'pedigree_data': self.pedigree_data, 'user_id': 'test_XXX'}
        self.client.force_authenticate(user=BwsTests.user)
        response = self.client.post(BwsTests.url, data, format='multipart',
                                    HTTP_ACCEPT="application/json")
        self.assertEqual(response.status_code, status.HTTP_408_REQUEST_TIMEOUT)
        content = json.loads(force_text(response.content))
        self.assertTrue('detail' in content)
        self.assertTrue('Request has timed out.' in content['detail'])
