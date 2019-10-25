""" Ovarian web-service testing.  """

from django.contrib.auth.models import User, Permission
from django.core.urlresolvers import reverse
from django.test import TestCase
from django.utils.encoding import force_text
from rest_framework import status
from rest_framework.authtoken.models import Token
from rest_framework.test import APIClient
import json
import os
from django.test.utils import override_settings


class OwsTests(TestCase):
    TEST_BASE_DIR = os.path.dirname(os.path.dirname(__file__))
    TEST_DATA_DIR = os.path.join(TEST_BASE_DIR, 'tests', 'data')

    @classmethod
    def setUpClass(cls):
        ''' Create a user and set up the test client. '''
        super(OwsTests, cls).setUpClass()
        cls.client = APIClient(enforce_csrf_checks=True)
        cls.user = User.objects.create_user('testuser', email='testuser@test.com',
                                            password='testing')
        # add user details
        # UserDetails.objects.create(user=cls.user, job_title=UserDetails.CGEN,
        #                            country='UK')
        # cls.user.user_permissions.add(Permission.objects.get(name='Can risk'))
        cls.user.save()
        cls.token = Token.objects.create(user=cls.user)
        cls.token.save()
        cls.url = reverse('ows')

    def setUp(self):
        self.pedigree_data = open(os.path.join(OwsTests.TEST_DATA_DIR, "canrisk_data_v1.txt"), "r")

    def tearDown(self):
        TestCase.tearDown(self)
        self.pedigree_data.close()

    def test_ows_output(self):
        ''' Test output of POSTing to the OWS using token authentication. '''
        data = {'mut_freq': 'UK', 'cancer_rates': 'UK',
                'pedigree_data': self.pedigree_data,
                'user_id': 'test_XXX'}
        OwsTests.client.credentials(HTTP_AUTHORIZATION='Token ' + OwsTests.token.key)
        response = OwsTests.client.post(OwsTests.url, data, format='multipart', HTTP_ACCEPT="application/json")
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        content = json.loads(force_text(response.content))
        self.assertTrue("mutation_frequency" in content)
        self.assertTrue("pedigree_result" in content)

        # check results returned for input pedigree
        pedigree_result = content["pedigree_result"][0]
        self.assertTrue("cancer_risks" in pedigree_result)
        self.assertGreater(len(pedigree_result["cancer_risks"]), 0)
        self.assertTrue("family_id" in pedigree_result)

    def test_multi_pedigree_ows(self):
        ''' Test POSTing multiple pedigrees to the OWS. '''
        multi_pedigree_data = open(os.path.join(OwsTests.TEST_DATA_DIR, "multi_canrisk_data_v1.txt"), "r")
        data = {'mut_freq': 'UK', 'cancer_rates': 'UK',
                'pedigree_data': multi_pedigree_data,
                'user_id': 'test_XXX'}
        OwsTests.client.credentials(HTTP_AUTHORIZATION='Token ' + OwsTests.token.key)
        response = OwsTests.client.post(OwsTests.url, data, format='multipart',
                                        HTTP_ACCEPT="application/json")
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        content = json.loads(force_text(response.content))
        self.assertEqual(len(content['pedigree_result']), 3, "three results")
        family_ids = ["XXX0", "XXX1", "XXX2"]
        for res in content['pedigree_result']:
            self.assertTrue(res['family_id'] in family_ids)

    def test_ows_warnings(self):
        ''' Test warning when proband has already had ovarian cancer and no risks are reported. '''
        # change proband to have had OC
        pd = self.pedigree_data.read().replace('1967\t0\t0\t0', '1967\t0\t0\t51')
        data = {'mut_freq': 'UK', 'cancer_rates': 'UK',
                'pedigree_data': pd, 'user_id': 'test_XXX'}
        OwsTests.client.force_authenticate(user=OwsTests.user)
        response = OwsTests.client.post(OwsTests.url, data, format='multipart',
                                        HTTP_ACCEPT="application/json")
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        content = json.loads(force_text(response.content))
        self.assertTrue('cancer_risks not provided' in content['warnings'])

    @override_settings(FORTRAN_TIMEOUT=0.01)
    def test_ows_timeout(self):
        ''' Test a timeout error is reported by the web-service. '''
        data = {'mut_freq': 'UK', 'cancer_rates': 'UK',
                'pedigree_data': self.pedigree_data, 'user_id': 'test_XXX'}
        OwsTests.client.force_authenticate(user=OwsTests.user)
        response = OwsTests.client.post(OwsTests.url, data, format='multipart',
                                        HTTP_ACCEPT="application/json")
        self.assertEqual(response.status_code, status.HTTP_408_REQUEST_TIMEOUT)
        content = json.loads(force_text(response.content))
        self.assertTrue('detail' in content)
        self.assertTrue('Request has timed out.' in content['detail'])


class OwsTestsPRS(TestCase):

    @classmethod
    def setUpClass(cls):
        ''' Create a user and set up the test client. '''
        super(OwsTestsPRS, cls).setUpClass()
        cls.client = APIClient(enforce_csrf_checks=True)
        cls.user = User.objects.create_user('testuser', email='testuser@test.com',
                                            password='testing')
        # add user details
        # UserDetails.objects.create(user=cls.user, job_title=UserDetails.CGEN,
        #                            country='UK')
        cls.user.user_permissions.add(Permission.objects.get(name='Can risk'))
        cls.user.save()
        cls.token = Token.objects.create(user=cls.user)
        cls.token.save()
        cls.url = reverse('ows')

    def setUp(self):
        self.pedigree_data = open(os.path.join(OwsTests.TEST_DATA_DIR, "canrisk_data_v1.txt"), "r")

    def test_prs_in_canrisk_file(self):
        '''
        Test ovarian cancer PRS parameters defined in the header of CanRisk formatted file.
        Calculate the cancer risk with and without a PRS and ensure that they are different.
        '''
        data = {'mut_freq': 'UK', 'cancer_rates': 'UK', 'pedigree_data': self.pedigree_data, 'user_id': 'test_XXX'}
        OwsTestsPRS.client.credentials(HTTP_AUTHORIZATION='Token ' + OwsTestsPRS.token.key)
        response = OwsTestsPRS.client.post(OwsTestsPRS.url, data, format='multipart', HTTP_ACCEPT="application/json")
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        orisk1 = json.loads(force_text(response.content))

        ped = open(os.path.join(OwsTests.TEST_DATA_DIR, "canrisk_data_v1.txt"), "r")
        pd = ped.read().replace('##CanRisk 1.0', '##CanRisk 1.0\n##PRS_OC=alpha=0.45,zscore=0.982')
        data = {'mut_freq': 'UK', 'cancer_rates': 'UK', 'pedigree_data': pd, 'user_id': 'test_XXX'}
        response = OwsTestsPRS.client.post(OwsTestsPRS.url, data, format='multipart', HTTP_ACCEPT="application/json")
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        orisk2 = json.loads(force_text(response.content))

        self.assertNotEqual(self.get_percent(orisk1, 80), self.get_percent(orisk2, 80),
                            "ovarian cancer at 80 different values")
        ped.close()

        # test with DEPRECATED beta instead of zscore
        ped = open(os.path.join(OwsTests.TEST_DATA_DIR, "canrisk_data_v1.txt"), "r")
        pd = ped.read().replace('##CanRisk 1.0', '##CanRisk 1.0\n##PRS_OC=alpha=0.45,beta=0.982')
        data = {'mut_freq': 'UK', 'cancer_rates': 'UK', 'pedigree_data': pd, 'user_id': 'test_XXX'}
        response = OwsTestsPRS.client.post(OwsTestsPRS.url, data, format='multipart', HTTP_ACCEPT="application/json")
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        orisk3 = json.loads(force_text(response.content))

        self.assertEqual(self.get_percent(orisk3, 80), self.get_percent(orisk2, 80),
                         "ovarian cancer at 80 different values")

    def get_percent(self, content, age):
        ''' Utility to return cancer percentage given the response content and an age. '''
        crisks = content["pedigree_result"][0]["cancer_risks"]
        for or2 in crisks:
            if or2['age'] == age:
                return or2['ovarian cancer risk']['percent']
