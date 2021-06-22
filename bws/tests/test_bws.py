""" BOADICEA web-service testing.  """
import os

from django.conf import settings
from django.contrib.auth.models import User
from django.urls import reverse
from django.test import TestCase
from rest_framework import status
from rest_framework.authtoken.models import Token
from rest_framework.test import APIClient
from django.utils.encoding import force_text
import json
from django.test.utils import override_settings


class BwsMixin(TestCase):
    TEST_BASE_DIR = os.path.dirname(os.path.dirname(__file__))
    TEST_DATA_DIR = os.path.join(TEST_BASE_DIR, 'tests', 'data')

    @classmethod
    def setUpClass(cls):
        ''' Create a user, token and url. '''
        super(BwsMixin, cls).setUpClass()
        cls.user = User.objects.create_user('testuser', email='testuser@test.com', password='testing')
        cls.user.save()
        cls.token = Token.objects.create(user=cls.user)
        cls.token.save()
        cls.url = reverse('bws')
        cls.client = APIClient(enforce_csrf_checks=True)
        cls.client.credentials(HTTP_AUTHORIZATION='Token ' + cls.token.key)


class MutFreqTests(BwsMixin):

    def test_ashkn_mut_freq(self):
        '''
        Test POSTing CanRisk file with multiple families with and without Ashkenazi Jewish ancestry.
        Check mutation frequencies are correctly set in each case.
        '''
        pedigree_data = open(os.path.join(BwsTests.TEST_DATA_DIR, "canrisk_multi4xAJ.txt"), "r")
        data = {'mut_freq': 'UK', 'cancer_rates': 'UK', 'pedigree_data': pedigree_data, 'user_id': 'test_XXX'}
        response = MutFreqTests.client.post(MutFreqTests.url, data, format='multipart', HTTP_ACCEPT="application/json")
        pedigree_data.close()
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        content = json.loads(force_text(response.content))
        self.assertTrue("mutation_frequency" in content)
        self.assertTrue("pedigree_result" in content)
        self.assertEqual(len(content['pedigree_result']), 4, "four results")

        family_ids = ["XXX2", "XX_AJ", "XX_NAJ", "XXX3"]
        for res in content['pedigree_result']:
            self.assertTrue(res['family_id'] in family_ids)

            if res['family_id'] == "XX_AJ" or res['family_id'] == "XX_NAJ":
                for r in res['cancer_risks']:
                    if r['age'] == 80:
                        c80 = r['breast cancer risk']['decimal']
                if res['family_id'] == "XX_AJ":
                    self.assertTrue("Ashkenazi" in res['mutation_frequency'])
                    c80_aj1 = c80
                elif res['family_id'] == "XX_NAJ":
                    self.assertTrue("UK" in res['mutation_frequency'])
                    c80_aj2 = c80

        self.assertGreater(c80_aj1, c80_aj2)
        self.assertTrue('mutation frequencies set to Ashkenazi Jewish population values for family (XX_AJ) as a ' +
                        'family member has Ashkenazi Jewish status.' in content['warnings'])


class BwsTests(BwsMixin):

    def setUp(self):
        ''' Set up test client and pedigree data. '''
        self.pedigree_data = open(os.path.join(BwsTests.TEST_DATA_DIR, "pedigree_data.txt"), "r")

    def tearDown(self):
        TestCase.tearDown(self)
        self.pedigree_data.close()

    def test_token_auth_bws(self):
        ''' Test POSTing to the BWS using token authentication. '''
        data = {'mut_freq': 'UK', 'cancer_rates': 'UK', 'pedigree_data': self.pedigree_data, 'user_id': 'test_XXX'}
        response = BwsTests.client.post(BwsTests.url, data, format='multipart', HTTP_ACCEPT="application/json")
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        content = json.loads(force_text(response.content))
        self.assertTrue("mutation_frequency" in content)
        self.assertTrue("pedigree_result" in content)
        self.assertTrue("family_id" in content["pedigree_result"][0])

    def test_multi_pedigree_bws(self):
        ''' Test POSTing multiple pedigrees to the BWS. '''
        multi_pedigree_data = open(os.path.join(BwsTests.TEST_DATA_DIR, "multi_pedigree_data.txt"), "r")
        data = {'mut_freq': 'UK', 'cancer_rates': 'UK', 'pedigree_data': multi_pedigree_data, 'user_id': 'test_XXX'}
        response = BwsTests.client.post(BwsTests.url, data, format='multipart', HTTP_ACCEPT="application/json")
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        content = json.loads(force_text(response.content))
        self.assertEqual(len(content['pedigree_result']), 2, "two results")
        family_ids = ["XXX0", "XXX1"]
        for res in content['pedigree_result']:
            self.assertTrue(res['family_id'] in family_ids)
        multi_pedigree_data.close()

    def test_canrisk_format_bws(self):
        ''' Test POSTing canrisk format pedigree to the BWS. '''
        canrisk_data = open(os.path.join(BwsTests.TEST_DATA_DIR, "canrisk_data_v1.txt"), "r")
        data = {'mut_freq': 'UK', 'cancer_rates': 'UK', 'pedigree_data': canrisk_data, 'user_id': 'test_XXX'}
        response = BwsTests.client.post(BwsTests.url, data, format='multipart', HTTP_ACCEPT="application/json")
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        content = json.loads(force_text(response.content))
        self.assertTrue("pedigree_result" in content)
        genes = settings.BC_MODEL['GENES']
        for g in genes:
            self.assertTrue(g in content['mutation_frequency']['UK'])
        canrisk_data.close()

    def test_token_auth_err(self):
        ''' Test POSTing to the BWS using token authentication. '''
        data = {'mut_freq': 'UK', 'cancer_rates': 'UK', 'pedigree_data': self.pedigree_data}
        client = APIClient(enforce_csrf_checks=True)
        client.credentials(HTTP_AUTHORIZATION='Token xxxxxxxxxx')
        response = client.post(BwsTests.url, data, format='multipart', HTTP_ACCEPT="application/json")
        self.assertNotEqual(response.status_code, status.HTTP_200_OK)
        content = json.loads(force_text(response.content))
        self.assertEqual(content['detail'], 'Invalid token.')

    def test_force_auth_bws(self):
        ''' Test POSTing to the BWS bypassing authentication. '''
        data = {'mut_freq': 'UK', 'cancer_rates': 'UK', 'pedigree_data': self.pedigree_data, 'user_id': 'test_XXX'}
        client = APIClient(enforce_csrf_checks=True)
        client.force_authenticate(user=BwsTests.user)
        response = client.post(BwsTests.url, data, format='multipart')
        self.assertEqual(response.status_code, status.HTTP_200_OK)

    # def test_custom_mutation_frequency(self):
        # ''' Test POSTing custom mutation frequencies. '''
        # genes = settings.BC_MODEL['GENES']
        # data = {g.lower() + '_mut_frequency': 0.00085 for g in genes}
        # data.update({'mut_freq': 'Custom', 'cancer_rates': 'UK',
        #              'pedigree_data': self.pedigree_data, 'user_id': 'test_XXX'})

        # response = BwsTests.client.post(BwsTests.url, data, format='multipart', HTTP_ACCEPT="application/json")
        # self.assertEqual(response.status_code, status.HTTP_200_OK)
        # content = json.loads(force_text(response.content))
        #
        # for g, mf in content['mutation_frequency']['Custom'].items():
        #     self.assertTrue(g in genes)
        #     self.assertEqual(mf, data[g.lower() + '_mut_frequency'])

    # def test_custom_mutation_frequency_errs(self):
        # ''' Test POSTing custom mutation frequencies with errors. '''
        # genes = settings.BC_MODEL['GENES']
        # data = {g.lower() + '_mut_frequency':
        #         (settings.MAX_MUTATION_FREQ + 0.1) if idx % 2 == 0 else (settings.MAX_MUTATION_FREQ - 0.1)
        #         for idx, g in enumerate(genes)}
        # data.update({'mut_freq': 'Custom', 'cancer_rates': 'UK',
        #              'pedigree_data': self.pedigree_data, 'user_id': 'test_XXX'})

        # response = BwsTests.client.post(BwsTests.url, data, format='multipart', HTTP_ACCEPT="application/json")
        # self.assertEqual(response.status_code, status.HTTP_400_BAD_REQUEST)
        # content = json.loads(force_text(response.content))
        # self.assertEqual(len(content.keys()), len(genes))
        # for k in content.keys():
        #     self.assertTrue(k.split("_")[0].upper() in genes)

    def test_missing_fields(self):
        ''' Test POSTing with missing required fields. '''
        data = {'mut_freq': 'UK'}
        response = BwsTests.client.post(BwsTests.url, data, format='multipart', HTTP_ACCEPT="application/json")
        self.assertEqual(response.status_code, status.HTTP_400_BAD_REQUEST)
        content = json.loads(force_text(response.content))
        self.assertEqual(content['user_id'][0], 'This field is required.')
        self.assertEqual(content['cancer_rates'][0], 'This field is required.')
        self.assertEqual(content['pedigree_data'][0], 'This field is required.')

    def test_bws_errors(self):
        ''' Test an error is reported by the web-service for an invalid year of birth. '''
        # force an error changing to an invalid year of birth
        ped = open(os.path.join(BwsTests.TEST_DATA_DIR, "pedigree_data.txt"), "r")
        pd = ped.read().replace('1963', '1600')
        data = {'mut_freq': 'UK', 'cancer_rates': 'UK', 'pedigree_data': pd, 'user_id': 'test_XXX'}
        response = BwsTests.client.post(BwsTests.url, data, format='multipart', HTTP_ACCEPT="application/json")
        self.assertEqual(response.status_code, status.HTTP_400_BAD_REQUEST)
        content = json.loads(force_text(response.content))
        ped.close()
        self.assertTrue('Person Error' in content)
        self.assertTrue('year of birth' in content['Person Error'])

    @override_settings(FORTRAN_TIMEOUT=0.05)
    def test_bws_timeout(self):
        ''' Test a timeout error is reported by the web-service. '''
        data = {'mut_freq': 'UK', 'cancer_rates': 'UK', 'pedigree_data': self.pedigree_data, 'user_id': 'test_XXX'}
        response = BwsTests.client.post(BwsTests.url, data, format='multipart', HTTP_ACCEPT="application/json")
        self.assertEqual(response.status_code, status.HTTP_408_REQUEST_TIMEOUT)
        content = json.loads(force_text(response.content))
        self.assertTrue('detail' in content)
        self.assertTrue('Request has timed out.' in content['detail'])


class TenYrTests(BwsMixin):

    def test_multi_tenyr(self):
        ''' Test POSTing multiple ages to the 10 year web service. '''

        tenyr_ages = "[30,40,45]"
        canrisk_data = open(os.path.join(TenYrTests.TEST_DATA_DIR, "canrisk_data_v1.txt"), "r")
        data = {'mut_freq': 'UK', 'cancer_rates': 'UK', 'pedigree_data': canrisk_data,
                'tenyr_ages': tenyr_ages, 'user_id': 'test_XXX'}

        response = TenYrTests.client.post(reverse('bcten'), data, format='multipart', HTTP_ACCEPT="application/json")
        canrisk_data.close()
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        content = json.loads(force_text(response.content))
        ten_yr_risks = content["pedigree_result"][0]["ten_yr_cancer_risk"]
        self.assertEqual(len(json.loads(tenyr_ages)), len(ten_yr_risks), "No. ranges equals no. 10 year risks")

    def test_tenyr_40(self):
        ''' Test 10 year risk from 40 same as BWS. '''

        bcten_url = reverse('bcten')
        bws_url = reverse('bws')
        dfs = [os.path.join(TenYrTests.TEST_DATA_DIR, "canrisk_data_v1.txt"),
               os.path.join(TenYrTests.TEST_DATA_DIR, "canrisk_data_v2.txt")]

        for df in dfs:
            # 1. 10 yr risk web-service for 40-50
            tenyr_ages = "[40]"
            canrisk_data = open(df, "r")
            data = {'mut_freq': 'UK', 'cancer_rates': 'UK', 'pedigree_data': canrisk_data,
                    'tenyr_ages': tenyr_ages, 'user_id': 'test_XXX'}
            response = TenYrTests.client.post(bcten_url, data, format='multipart', HTTP_ACCEPT="application/json")
            canrisk_data.close()
            self.assertEqual(response.status_code, status.HTTP_200_OK)
            content = json.loads(force_text(response.content))
            ten_yr_cancer_risk1 = content["pedigree_result"][0]["ten_yr_cancer_risk"][0]

            # 2. use BWS to get 10 yr risk from 40 for comparison
            canrisk_data = open(df, "r")
            data = {'mut_freq': 'UK', 'cancer_rates': 'UK', 'pedigree_data': canrisk_data, 'user_id': 'test_XXX'}
            response = TenYrTests.client.post(bws_url, data, format='multipart', HTTP_ACCEPT="application/json")
            canrisk_data.close()
            self.assertEqual(response.status_code, status.HTTP_200_OK)
            content = json.loads(force_text(response.content))
            ten_yr_cancer_risk2 = content["pedigree_result"][0]['ten_yr_cancer_risk'][0]
            self.assertDictEqual(ten_yr_cancer_risk1, ten_yr_cancer_risk2, "Compare 10yr cancer risk from 40")
