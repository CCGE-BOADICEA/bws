"""
Prostate web-service testing.

© 2024 University of Cambridge
SPDX-FileCopyrightText: 2024 University of Cambridge
SPDX-License-Identifier: GPL-3.0-or-later
"""

from django.conf import settings
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
import unittest


class PwsTests(TestCase):
    TEST_BASE_DIR = os.path.dirname(os.path.dirname(__file__))
    TEST_DATA_DIR = os.path.join(TEST_BASE_DIR, 'tests', 'data')

    @classmethod
    def setUpClass(cls):
        ''' Create a user and set up the test client. '''
        super(PwsTests, cls).setUpClass()
        if not settings.PROSTATE_CANCER:
            return

        cls.client = APIClient(enforce_csrf_checks=True)
        cls.user = User.objects.create_user('testuser', email='testuser@test.com',
                                            password='testing')
        cls.user.user_permissions.add(Permission.objects.get(name='Commercial PC webservices'))
        cls.user.save()
        cls.token = Token.objects.create(user=cls.user)
        cls.token.save()
        cls.url = reverse('pws')

    def setUp(self):
        self.pedigree_data = open(os.path.join(PwsTests.TEST_DATA_DIR, "male.canrisk3"), "r")

    def tearDown(self):
        TestCase.tearDown(self)
        self.pedigree_data.close()

    @unittest.skipIf(not settings.PROSTATE_CANCER, "prostate cancer model not used")
    def test_pws_output(self):
        ''' Test output of POSTing to the PWS using token authentication. '''
        data = {'mut_freq': 'UK', 'cancer_rates': 'UK', 'pedigree_data': self.pedigree_data,
                'user_id': 'test_XXX', 'prs': json.dumps({'alpha': 0.45, 'zscore': 1.652})}
        PwsTests.client.credentials(HTTP_AUTHORIZATION='Token ' + PwsTests.token.key)
        response = PwsTests.client.post(PwsTests.url, data, format='multipart', HTTP_ACCEPT="application/json")
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        content = json.loads(force_str(response.content))
        self.assertTrue("mutation_frequency" in content)
        self.assertTrue("pedigree_result" in content)

        # check results returned for input pedigree
        pedigree_result = content["pedigree_result"][0]
        self.assertTrue("cancer_risks" in pedigree_result)
        self.assertGreater(len(pedigree_result["cancer_risks"]), 0)
        self.assertTrue("family_id" in pedigree_result)

    @unittest.skipIf(not settings.PROSTATE_CANCER, "prostate cancer model not used")
    def test_pws_warnings(self):
        ''' Test warning when proband has already had prostate cancer and no risks are reported. '''
        # change proband to have had OC
        pd = self.pedigree_data.read().replace('1962\t0\t0\t0\t0', '1962\t0\t0\t0\t51')
        data = {'mut_freq': 'UK', 'cancer_rates': 'UK', 'pedigree_data': pd, 'user_id': 'test_XXX'}
        PwsTests.client.force_authenticate(user=PwsTests.user)
        response = PwsTests.client.post(PwsTests.url, data, format='multipart', HTTP_ACCEPT="application/json")
        self.assertEqual(response.status_code, status.HTTP_400_BAD_REQUEST)
        content = json.loads(force_str(response.content))
        self.assertTrue('PROBAND HAS ALREADY HAD A CANCER' in content['Model Error'])

    @unittest.skipIf(not settings.PROSTATE_CANCER, "prostate cancer model not used")
    @override_settings(FORTRAN_TIMEOUT=0.01)
    def test_pws_timeout(self):
        ''' Test a timeout error is reported by the web-service. '''
        data = {'mut_freq': 'UK', 'cancer_rates': 'UK',
                'pedigree_data': self.pedigree_data, 'user_id': 'test_XXX'}
        PwsTests.client.force_authenticate(user=PwsTests.user)
        response = PwsTests.client.post(PwsTests.url, data, format='multipart', HTTP_ACCEPT="application/json")
        self.assertEqual(response.status_code, status.HTTP_408_REQUEST_TIMEOUT)
        content = json.loads(force_str(response.content))
        self.assertTrue('detail' in content)
        self.assertTrue('Request has timed out.' in content['detail'])


class PwsTestsPRS(TestCase):

    @classmethod
    def setUpClass(cls):
        ''' Create a user and set up the test client. '''
        super(PwsTestsPRS, cls).setUpClass()
        if not settings.PROSTATE_CANCER:
            return

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
        cls.url = reverse('pws')

    def setUp(self):
        self.pedigree_data = open(os.path.join(PwsTests.TEST_DATA_DIR, "male.canrisk3"), "r")

    @unittest.skipIf(not settings.PROSTATE_CANCER, "prostate cancer model not used")
    def test_prs_in_canrisk_file(self):
        '''
        Test prostate cancer PRS parameters defined in the header of CanRisk formatted file.
        Calculate the cancer risk with and without a PRS and ensure that they are different.
        '''
        data = {'mut_freq': 'UK', 'cancer_rates': 'UK', 'pedigree_data': self.pedigree_data, 'user_id': 'test_XXX'}
        PwsTestsPRS.client.credentials(HTTP_AUTHORIZATION='Token ' + PwsTestsPRS.token.key)
        response = PwsTestsPRS.client.post(PwsTestsPRS.url, data, format='multipart', HTTP_ACCEPT="application/json")
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        prisk1 = json.loads(force_str(response.content))

        ped = open(os.path.join(PwsTests.TEST_DATA_DIR, "male.canrisk3"), "r")
        pd = ped.read().replace('##CanRisk 3', '##CanRisk 4\n##PRS_PC=alpha=0.45,zscore=0.982')
        ped.close()
        pd = pd.replace("\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0",         # add extra column for HOXB13
                        "\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0")
        data = {'mut_freq': 'UK', 'cancer_rates': 'UK', 'pedigree_data': pd, 'user_id': 'test_XXX'}

        response = PwsTestsPRS.client.post(PwsTestsPRS.url, data, format='multipart', HTTP_ACCEPT="application/json")
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        prisk2 = json.loads(force_str(response.content))

        self.assertNotEqual(self.get_percent(prisk1, 80), self.get_percent(prisk2, 80),
                            "prostate cancer at 80 different values")

    def get_percent(self, content, age):
        ''' Utility to return cancer percentage given the response content and an age. '''
        crisks = content["pedigree_result"][0]["cancer_risks"]
        for r in crisks:
            if r['age'] == age:
                return r['prostate cancer risk']['percent']
