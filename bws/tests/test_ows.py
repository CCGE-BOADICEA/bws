"""
Ovarian web-service testing.

© 2023 University of Cambridge
SPDX-FileCopyrightText: 2023 University of Cambridge
SPDX-License-Identifier: GPL-3.0-or-later
"""

from django.contrib.auth.models import User, Permission
from django.test import TestCase
from django.test.utils import override_settings
from django.urls import reverse
from django.utils.encoding import force_str
from rest_framework import status
from rest_framework.authtoken.models import Token
from rest_framework.test import APIClient
import json
import os
from bws.pedigree_file import PedigreeFile


class OwsTests(TestCase):
    TEST_BASE_DIR = os.path.dirname(os.path.dirname(__file__))
    TEST_DATA_DIR = os.path.join(TEST_BASE_DIR, 'tests', 'data')

    @classmethod
    def setUpClass(cls):
        ''' Create a user and set up the test client. '''
        super(OwsTests, cls).setUpClass()
        cls.drf_client = APIClient(enforce_csrf_checks=True)
        cls.user = User.objects.create_user('testuser', email='testuser@test.com',
                                            password='testing')
        # add user details
        # UserDetails.objects.create(user=cls.user, job_title=UserDetails.CGEN,
        #                            country='UK')
        cls.user.user_permissions.add(Permission.objects.get(name='Commercial OC webservices'))
        cls.user.save()
        cls.token = Token.objects.create(user=cls.user)
        cls.token.save()
        cls.url = reverse('ows')

    def setUp(self):
        self.pedigree_data = open(os.path.join(OwsTests.TEST_DATA_DIR, "d0.canrisk"), "r")

    def tearDown(self):
        TestCase.tearDown(self)
        self.pedigree_data.close()

    def test_ows_output(self):
        ''' Test output of POSTing to the OWS using token authentication. '''
        data = {'mut_freq': 'UK', 'cancer_rates': 'UK',
                'pedigree_data': self.pedigree_data,
                'user_id': 'test_XXX', 'prs': json.dumps({'alpha': 0.45, 'zscore': 1.652})}
        OwsTests.drf_client.credentials(HTTP_AUTHORIZATION='Token ' + OwsTests.token.key)
        response = OwsTests.drf_client.post(OwsTests.url, data, format='multipart', HTTP_ACCEPT="application/json")
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        content = json.loads(force_str(response.content))
        self.assertTrue("mutation_frequency" in content)
        self.assertTrue("pedigree_result" in content)

        # check results returned for input pedigree
        pedigree_result = content["pedigree_result"][0]
        self.assertTrue("lifetime_cancer_risk" in pedigree_result)
        self.assertTrue("baseline_lifetime_cancer_risk" in pedigree_result)
        self.assertTrue("baseline_cancer_risks" in pedigree_result)
        self.assertTrue("cancer_risks" in pedigree_result)
        self.assertGreater(len(pedigree_result["cancer_risks"]), 0)
        self.assertTrue("family_id" in pedigree_result)

    def test_multi_pedigree_ows(self):
        ''' Test POSTing multiple pedigrees to the OWS. '''
        multi_pedigree_data = open(os.path.join(OwsTests.TEST_DATA_DIR, "multi", "d2.canrisk"), "r")
        data = {'mut_freq': 'UK', 'cancer_rates': 'UK',
                'pedigree_data': multi_pedigree_data,
                'user_id': 'test_XXX'}
        OwsTests.drf_client.credentials(HTTP_AUTHORIZATION='Token ' + OwsTests.token.key)
        response = OwsTests.drf_client.post(OwsTests.url, data, format='multipart', HTTP_ACCEPT="application/json")
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        content = json.loads(force_str(response.content))
        self.assertEqual(len(content['pedigree_result']), 3, "three results")
        multi_pedigree_data.close()

        #
        with open(os.path.join(OwsTests.TEST_DATA_DIR, "multi", "d2.canrisk"), "r") as f:
            multi_pedigree_data = f.read()
        f.close()
        pf = PedigreeFile(multi_pedigree_data)
        family_ids = []
        pedigrees = {}
        for p in pf.pedigrees:
            family_ids.append(p.famid)
            pedigrees[p.famid] = p

        for res in content['pedigree_result']:
            self.assertTrue(res['family_id'] in family_ids)
            p = pedigrees[res['family_id']]
            t = p.get_target()
            if not t.cancers.is_cancer_diagnosed(): # if no cancer diagnosis lifetime cancer risks present
                self.assertTrue("lifetime_cancer_risk" in res)
                self.assertTrue("baseline_lifetime_cancer_risk" in res)
            
            self.assertTrue("baseline_cancer_risks" in res)
            self.assertTrue("cancer_risks" in res)

    def test_ows_bwa(self):
        ''' Test web-service takes BWA file as input. '''
        pedigree_data = open(os.path.join(OwsTests.TEST_DATA_DIR, "d3.bwa"), "r")
        data = {'mut_freq': 'UK', 'cancer_rates': 'UK',
                'pedigree_data': pedigree_data,
                'user_id': 'test_XXX'}
        OwsTests.drf_client.credentials(HTTP_AUTHORIZATION='Token ' + OwsTests.token.key)
        response = OwsTests.drf_client.post(OwsTests.url, data, format='multipart', HTTP_ACCEPT="application/json")
        self.assertEqual(response.status_code, status.HTTP_200_OK)

    def test_ows_warnings(self):
        ''' Test warning when proband has already had ovarian cancer and no risks are reported. '''
        # change proband to have had OC
        pd = self.pedigree_data.read().replace('1967\t0\t0\t0', '1967\t0\t0\t51')
        data = {'mut_freq': 'UK', 'cancer_rates': 'UK',
                'pedigree_data': pd, 'user_id': 'test_XXX'}
        OwsTests.drf_client.force_authenticate(user=OwsTests.user)
        response = OwsTests.drf_client.post(OwsTests.url, data, format='multipart', HTTP_ACCEPT="application/json")
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        content = json.loads(force_str(response.content))
        self.assertTrue('cancer_risks not provided' in content['warnings'])

    @override_settings(FORTRAN_TIMEOUT=0.001)
    def test_ows_timeout(self):
        ''' Test a timeout error is reported by the web-service. '''
        data = {'mut_freq': 'UK', 'cancer_rates': 'UK',
                'pedigree_data': self.pedigree_data, 'user_id': 'test_XXX'}
        OwsTests.drf_client.force_authenticate(user=OwsTests.user)
        response = OwsTests.drf_client.post(OwsTests.url, data, format='multipart', HTTP_ACCEPT="application/json")
        self.assertEqual(response.status_code, status.HTTP_408_REQUEST_TIMEOUT)
        content = json.loads(force_str(response.content))
        self.assertTrue('detail' in content)
        self.assertTrue('Request has timed out.' in content['detail'])


class OwsTestsPRS(TestCase):

    @classmethod
    def setUpClass(cls):
        ''' Create a user and set up the test client. '''
        super(OwsTestsPRS, cls).setUpClass()
        cls.drf_client = APIClient(enforce_csrf_checks=True)
        cls.user = User.objects.create_user('testuser', email='testuser@test.com',
                                            password='testing')
        # add user details
        # UserDetails.objects.create(user=cls.user, job_title=UserDetails.CGEN,
        #                            country='UK')
        cls.user.user_permissions.add(Permission.objects.get(name='Commercial OC webservices'))
        cls.user.save()
        cls.token = Token.objects.create(user=cls.user)
        cls.token.save()
        cls.url = reverse('ows')

    def setUp(self):
        self.pedigree_data = open(os.path.join(OwsTests.TEST_DATA_DIR, "d0.canrisk"), "r")

    def test_prs_in_canrisk_file(self):
        '''
        Test ovarian cancer PRS parameters defined in the header of CanRisk formatted file.
        Calculate the cancer risk with and without a PRS and ensure that they are different.
        '''
        data = {'mut_freq': 'UK', 'cancer_rates': 'UK', 'pedigree_data': self.pedigree_data, 'user_id': 'test_XXX'}
        OwsTestsPRS.drf_client.credentials(HTTP_AUTHORIZATION='Token ' + OwsTestsPRS.token.key)
        response = OwsTestsPRS.drf_client.post(OwsTestsPRS.url, data, format='multipart', HTTP_ACCEPT="application/json")
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        orisk1 = json.loads(force_str(response.content))

        ped = open(os.path.join(OwsTests.TEST_DATA_DIR, "d0.canrisk"), "r")
        pd = ped.read().replace('##CanRisk 1.0', '##CanRisk 1.0\n##PRS_OC=alpha=0.45,zscore=0.982')
        data = {'mut_freq': 'UK', 'cancer_rates': 'UK', 'pedigree_data': pd, 'user_id': 'test_XXX'}
        response = OwsTestsPRS.drf_client.post(OwsTestsPRS.url, data, format='multipart', HTTP_ACCEPT="application/json")
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        orisk2 = json.loads(force_str(response.content))

        self.assertNotEqual(self.get_percent(orisk1, 80), self.get_percent(orisk2, 80),
                            "ovarian cancer at 80 different values")
        ped.close()

        # test with DEPRECATED beta instead of zscore
        ped = open(os.path.join(OwsTests.TEST_DATA_DIR, "d0.canrisk"), "r")
        pd = ped.read().replace('##CanRisk 1.0', '##CanRisk 1.0\n##PRS_OC=alpha=0.45,beta=0.982')
        data = {'mut_freq': 'UK', 'cancer_rates': 'UK', 'pedigree_data': pd, 'user_id': 'test_XXX'}
        response = OwsTestsPRS.drf_client.post(OwsTestsPRS.url, data, format='multipart', HTTP_ACCEPT="application/json")
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        orisk3 = json.loads(force_str(response.content))

        self.assertEqual(self.get_percent(orisk3, 80), self.get_percent(orisk2, 80),
                         "ovarian cancer at 80 different values")

    def get_percent(self, content, age):
        ''' Utility to return cancer percentage given the response content and an age. '''
        crisks = content["pedigree_result"][0]["cancer_risks"]
        for or2 in crisks:
            if or2['age'] == age:
                return or2['ovarian cancer risk']['percent']
