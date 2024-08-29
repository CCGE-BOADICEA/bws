"""
BOADICEA web-service testing.

© 2023 University of Cambridge
SPDX-FileCopyrightText: 2023 University of Cambridge
SPDX-License-Identifier: GPL-3.0-or-later
"""

from bws.calc.calcs import Predictions
from bws.cancer import CanRiskGeneticTests
from bws.exceptions import ModelError
from bws.pedigree import CanRiskPedigree, Female
from datetime import date
from django.conf import settings
from django.contrib.auth.models import User, Permission
from django.test import TestCase
from django.test.utils import override_settings
from django.urls import reverse
from django.utils.encoding import force_str
from rest_framework import status
from rest_framework.authtoken.models import Token
from rest_framework.exceptions import ValidationError
from rest_framework.test import APIClient
import json
import os


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
        cls.user.user_permissions.add(Permission.objects.get(name='Can risk'))
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
        pedigree_data = open(os.path.join(BwsTests.TEST_DATA_DIR, "multi", "d3.4xAJ.canrisk2"), "r")
        data = {'mut_freq': 'UK', 'cancer_rates': 'UK', 'pedigree_data': pedigree_data, 'user_id': 'test_XXX'}
        response = MutFreqTests.client.post(MutFreqTests.url, data, format='multipart', HTTP_ACCEPT="application/json")
        pedigree_data.close()
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        content = json.loads(force_str(response.content))
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
        self.pedigree_data = open(os.path.join(BwsTests.TEST_DATA_DIR, "d3.bwa"), "r")

    def tearDown(self):
        TestCase.tearDown(self)
        self.pedigree_data.close()

    def test_token_auth_bws(self):
        ''' Test POSTing to the BWS using token authentication. '''
        data = {'mut_freq': 'UK', 'cancer_rates': 'UK', 'pedigree_data': self.pedigree_data, 'user_id': 'test_XXX'}
        response = BwsTests.client.post(BwsTests.url, data, format='multipart', HTTP_ACCEPT="application/json")
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        content = json.loads(force_str(response.content))
        self.assertTrue("mutation_frequency" in content)
        self.assertTrue("pedigree_result" in content)
        self.assertTrue("family_id" in content["pedigree_result"][0])

    def test_multi_pedigree_bws(self):
        ''' Test POSTing multiple pedigrees to the BWS. '''
        multi_pedigree_data = open(os.path.join(BwsTests.TEST_DATA_DIR, "multi", "d1.bwa"), "r")
        data = {'mut_freq': 'UK', 'cancer_rates': 'UK', 'pedigree_data': multi_pedigree_data, 'user_id': 'test_XXX'}
        response = BwsTests.client.post(BwsTests.url, data, format='multipart', HTTP_ACCEPT="application/json")
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        content = json.loads(force_str(response.content))
        self.assertEqual(len(content['pedigree_result']), 2, "two results")
        family_ids = ["XXX0", "XXX1"]
        for res in content['pedigree_result']:
            self.assertTrue(res['family_id'] in family_ids)
        multi_pedigree_data.close()

    def test_canrisk_format_bws(self):
        ''' Test POSTing canrisk format pedigree to the BWS. '''
        canrisk_data = open(os.path.join(BwsTests.TEST_DATA_DIR, "d0.canrisk"), "r")
        data = {'mut_freq': 'UK', 'cancer_rates': 'UK', 'pedigree_data': canrisk_data, 'user_id': 'test_XXX'}
        response = BwsTests.client.post(BwsTests.url, data, format='multipart', HTTP_ACCEPT="application/json")
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        content = json.loads(force_str(response.content))
        self.assertTrue("pedigree_result" in content)
        genes = settings.BC_MODEL['GENES']
        for g in genes:
            self.assertTrue(g in content['mutation_frequency']['UK'])
        canrisk_data.close()

    def test_canrisk_v2_format_bws(self):
        ''' Test POSTing canrisk format pedigree to the BWS. '''
        canrisk_data = open(os.path.join(BwsTests.TEST_DATA_DIR, "d1.canrisk2"), "r")
        data = {'mut_freq': 'UK', 'cancer_rates': 'France', 'pedigree_data': canrisk_data, 'user_id': 'test_XXX'}
        response = BwsTests.client.post(BwsTests.url, data, format='multipart', HTTP_ACCEPT="application/json")
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        content = json.loads(force_str(response.content))
        self.assertTrue("pedigree_result" in content)
        genes = settings.BC_MODEL['GENES']
        for g in genes:
            self.assertTrue(g in content['mutation_frequency']['UK'])

        self.assertDictEqual(settings.BC_MODEL['GENETIC_TEST_SENSITIVITY'], content['mutation_sensitivity'])
        self.assertEqual(content['cancer_incidence_rates'], 'France')
        canrisk_data.close()

    def test_token_auth_err(self):
        ''' Test POSTing to the BWS using token authentication. '''
        data = {'mut_freq': 'UK', 'cancer_rates': 'UK', 'pedigree_data': self.pedigree_data}
        client = APIClient(enforce_csrf_checks=True)
        client.credentials(HTTP_AUTHORIZATION='Token xxxxxxxxxx')
        response = client.post(BwsTests.url, data, format='multipart', HTTP_ACCEPT="application/json")
        self.assertNotEqual(response.status_code, status.HTTP_200_OK)
        content = json.loads(force_str(response.content))
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
        # content = json.loads(force_str(response.content))
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
        # content = json.loads(force_str(response.content))
        # self.assertEqual(len(content.keys()), len(genes))
        # for k in content.keys():
        #     self.assertTrue(k.split("_")[0].upper() in genes)

    def test_missing_fields(self):
        ''' Test POSTing with missing required fields. '''
        data = {'mut_freq': 'UK'}
        response = BwsTests.client.post(BwsTests.url, data, format='multipart', HTTP_ACCEPT="application/json")
        self.assertEqual(response.status_code, status.HTTP_400_BAD_REQUEST)
        content = json.loads(force_str(response.content))
        self.assertEqual(content['user_id'][0], 'This field is required.')
        self.assertEqual(content['cancer_rates'][0], 'This field is required.')
        self.assertEqual(content['pedigree_data'][0], 'No file was submitted.')

    def test_bws_errors(self):
        ''' Test an error is reported by the web-service for an invalid year of birth. '''
        # force an error changing to an invalid year of birth
        ped = open(os.path.join(BwsTests.TEST_DATA_DIR, "d3.bwa"), "r")
        pd = ped.read().replace('1963', '1600')
        data = {'mut_freq': 'UK', 'cancer_rates': 'UK', 'pedigree_data': pd, 'user_id': 'test_XXX'}
        response = BwsTests.client.post(BwsTests.url, data, format='multipart', HTTP_ACCEPT="application/json")
        self.assertEqual(response.status_code, status.HTTP_400_BAD_REQUEST)
        content = json.loads(force_str(response.content))
        ped.close()
        self.assertTrue('Person Error' in content)
        self.assertTrue('year of birth' in content['Person Error'])

    def test_field_validate_errors(self):
        ''' Test error with superfluous fields included. '''
        data = {'mut_freq': 'UK', 'nosuchflag': '1234', 'nosuchflag2': 77}
        response = BwsTests.client.post(BwsTests.url, data, format='multipart', HTTP_ACCEPT="application/json")
        self.assertEqual(response.status_code, status.HTTP_400_BAD_REQUEST)
        content = json.loads(force_str(response.content))
        self.assertTrue('Input Field Error' in content)
        self.assertEqual('Extra input field(s) found: nosuchflag, nosuchflag2', content['Input Field Error'])

    @override_settings(FORTRAN_TIMEOUT=0.0005)
    def test_bws_timeout(self):
        ''' Test a timeout error is reported by the web-service. '''
        data = {'mut_freq': 'UK', 'cancer_rates': 'UK', 'pedigree_data': self.pedigree_data, 'user_id': 'test_XXX'}
        response = BwsTests.client.post(BwsTests.url, data, format='multipart', HTTP_ACCEPT="application/json")
        self.assertEqual(response.status_code, status.HTTP_408_REQUEST_TIMEOUT)
        content = json.loads(force_str(response.content))
        self.assertTrue('detail' in content)
        self.assertTrue('Request has timed out.' in content['detail'])

    def test_calcs_validation_err(self):
        """ Test invalid calculation type raise ValidationError. """
        target = Female("FAM1", "F0", "001", "002", "003", target="1", age="20",
                        yob=str(date.today().year-20), gtests=CanRiskGeneticTests.default_factory())
        pedigree = CanRiskPedigree(people=[target])
        with self.assertRaisesRegex(ValidationError, r"Unknown calculation requested: dribble"):
            Predictions(pedigree, model_settings=settings.OC_MODEL, calcs=['dribble'])

    def test_bws_model_err(self):
        ''' Test ModelError raised because of pedigree file with only one twin. '''
        target = Female("FAM1", "F0", "001", "002", "003", target="1", age="20", mztwin="1",
                        yob=str(date.today().year-20), gtests=CanRiskGeneticTests.default_factory())
        pedigree = CanRiskPedigree(people=[target])
        with self.assertRaisesRegex(ModelError, r"ERRORS IN THE PEDIGREE FILE"):
            Predictions(pedigree)

    def test_bws_model_deceased_no_risks(self):
        ''' Test deceased target produces mutation carrier probabilities and no risks. '''
        canrisk_data = open(os.path.join(BwsTests.TEST_DATA_DIR, "d5.dead.canrisk2"), "r")
        data = {'mut_freq': 'UK', 'cancer_rates': 'France', 'pedigree_data': canrisk_data, 'user_id': 'test_XXX'}
        response = BwsTests.client.post(BwsTests.url, data, format='multipart', HTTP_ACCEPT="application/json")
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        content = json.loads(force_str(response.content))
        self.assertTrue("pedigree_result" in content)
        self.assertFalse('cancer_risks' in content["pedigree_result"][0])
        self.assertFalse('lifetime_cancer_risk' in content["pedigree_result"][0])
        self.assertFalse('baseline_cancer_risks' in content["pedigree_result"][0])
        genes = settings.BC_MODEL['GENES']
        for g in genes:
            self.assertTrue(g in content['mutation_frequency']['UK'])


class CombineModelResultsTests(BwsMixin):

    def test_results_page(self):
        ''' Get results from breast and ovarian web-services and test calling combine results web-service. '''
        # 1. calculate breast cancer risks
        fn = os.path.join(BwsTests.TEST_DATA_DIR, "d4.AJ.canrisk2")
        canrisk_data = open(fn, "r")
        data = {'mut_freq': 'UK', 'cancer_rates': 'Spain', 'pedigree_data': canrisk_data, 'user_id': 'test_XXX'}
        bws_result = CombineModelResultsTests.client.post(reverse('bws'), data, format='multipart',
                                                          HTTP_ACCEPT="application/json")
        self.assertEqual(bws_result.status_code, status.HTTP_200_OK)

        # 2. calculate ovarian cancer risks
        canrisk_data = open(fn, "r")
        data = {'mut_freq': 'UK', 'cancer_rates': 'Spain', 'pedigree_data': canrisk_data, 'user_id': 'test_XXX'}
        ows_result = CombineModelResultsTests.client.post(reverse('ows'), data,
                                                          format='multipart', HTTP_ACCEPT="application/json")
        self.assertEqual(ows_result.status_code, status.HTTP_200_OK)

        # 3. combine results
        data = {"ows_result": json.dumps(ows_result.json(), separators=(',', ':')),
                "bws_result": json.dumps(bws_result.json(), separators=(',', ':'))}

        response = CombineModelResultsTests.client.post(reverse('combine'), data=data)
        self.assertEqual(response.status_code, status.HTTP_200_OK)
