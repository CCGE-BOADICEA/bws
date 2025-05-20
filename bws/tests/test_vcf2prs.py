"""
© 2023 University of Cambridge
SPDX-FileCopyrightText: 2023 University of Cambridge
SPDX-License-Identifier: GPL-3.0-or-later
"""

import json
import os
from pathlib import Path

from django.contrib.auth.models import User, Permission
from django.test.testcases import TestCase
from django.test.utils import override_settings
from django.urls import reverse
from django.utils.encoding import force_str
from rest_framework import status
from rest_framework.authtoken.models import Token
from rest_framework.test import APIClient
import vcf2prs
from vcf2prs.prs import Prs
from bws.vcf2prs_api import Vcf2PrsView


class Vcf2PrsWebServices(TestCase):
    ''' Test the VCF to PRS webservice '''

    @classmethod
    def setUpClass(cls):
        super(Vcf2PrsWebServices, cls).setUpClass()
        cls.client = APIClient(enforce_csrf_checks=True)
        cls.user = User.objects.create_user('testuser', email='testuser@test.com',
                                            password='testing')
        # add user details
        # UserDetails.objects.create(user=cls.user, job_title=UserDetails.CGEN,
        #                            country='UK')
        cls.user.save()
        cls.user.user_permissions.add(Permission.objects.get(name='Commercial BC webservices'))
        cls.token = Token.objects.create(user=cls.user)
        cls.token.save()
        cls.url = reverse('prs')
        cls.moduledir = Path(vcf2prs.__file__).parent.parent

    def setUp(self):
        ''' Create a user and set up the test client. '''
        self.vcf_file = os.path.join(Vcf2PrsWebServices.moduledir, "Examples/VCFfiles/sample_BCAC_313.vcf")
        self.vcf_data = open(self.vcf_file, "r")
        self.prs_reference_file = "BCAC_313_PRS.prs"
        self.prs_file_name = os.path.join(Vcf2PrsWebServices.moduledir, "PRSmodels_CanRisk", self.prs_reference_file)

    def test_prs(self):
        ''' Test POSTing to a vcf file to get a polygenic risk score. '''
        data = {'vcf_file': self.vcf_data, 'sample_name': 't0.5', 'bc_prs_reference_file': self.prs_reference_file}
        Vcf2PrsWebServices.client.credentials(HTTP_AUTHORIZATION='Token ' + Vcf2PrsWebServices.token.key)
        response = Vcf2PrsWebServices.client.post(Vcf2PrsWebServices.url, data, format='multipart',
                                                  HTTP_ACCEPT="application/json")
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertContains(response, "alpha")
        self.assertContains(response, "zscore")

    def test_prs_v_direct(self):
        ''' Test POSTing to a vcf file to get a polygenic risk score and
        compare with the direct call to calculate a PRS. '''
        data = {'vcf_file': self.vcf_data, 'sample_name': 't0.9', 'bc_prs_reference_file': self.prs_reference_file}
        Vcf2PrsWebServices.client.credentials(HTTP_AUTHORIZATION='Token ' + Vcf2PrsWebServices.token.key)
        response = Vcf2PrsWebServices.client.post(Vcf2PrsWebServices.url, data, format='multipart',
                                                  HTTP_ACCEPT="application/json")
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        content = json.loads(force_str(response.content))

        prs = Prs(prs_file=self.prs_file_name)
        prs.calculate_prs_from_vcf(self.vcf_file, 't0.9')
        prs.calculate_z_from_raw(prs.raw_PRS)

        #prs = Prs(prs_file=self.prs_file_name, geno_file=self.vcf_file, sample='t0.9')
        #zscore = prs.z_Score
        self.assertEqual(prs.z_Score[0], content['breast_cancer_prs']['zscore'], 'web-service and direct calculation')

    def test_prs_err(self):
        ''' Test POSTing to the 400 returned without the vcf specified. '''
        data = {'sample_name': 'SampleA'}
        Vcf2PrsWebServices.client.credentials(HTTP_AUTHORIZATION='Token ' + Vcf2PrsWebServices.token.key)
        response = Vcf2PrsWebServices.client.post(Vcf2PrsWebServices.url, data, format='multipart',
                                                  HTTP_ACCEPT="application/json")
        self.assertEqual(response.status_code, status.HTTP_400_BAD_REQUEST)
        self.assertEqual(response.data['vcf_file'][0], 'No file was submitted.')

    @override_settings(DATA_UPLOAD_MAX_MEMORY_SIZE=1)
    def test_prs_upload_limit(self):
        ''' Test POSTing to a vcf file that is too large. '''
        data = {'vcf_file': self.vcf_data, 'sample_name': 'SampleA'}
        Vcf2PrsWebServices.client.credentials(HTTP_AUTHORIZATION='Token ' + Vcf2PrsWebServices.token.key)
        response = Vcf2PrsWebServices.client.post(Vcf2PrsWebServices.url, data, format='multipart',
                                                  HTTP_ACCEPT="application/json")
        self.assertEqual(response.status_code, status.HTTP_400_BAD_REQUEST)
