"""
© 2026 University of Cambridge
SPDX-FileCopyrightText: 2026 University of Cambridge
SPDX-License-Identifier: GPL-3.0-or-later
"""

from django.test.testcases import TestCase
from django.urls import reverse
from rest_framework import status
from rest_framework.test import APIClient


class AltchaWebServices(TestCase):
    ''' Test the altcha webservice '''

    def test_challenge_view_returns_200(self):
        ''' Test fetching a new random challenge to be used by the ALTCHA widget. '''
        drf_client = APIClient(enforce_csrf_checks=True)
        url = reverse('altcha')
        response = drf_client.get(url)
        self.assertEqual(response.status_code, status.HTTP_200_OK)

    def test_challenge_view_returns_json(self):
        ''' Test challenge response returns JSON. '''
        drf_client = APIClient(enforce_csrf_checks=True)
        url = reverse('altcha')
        response = drf_client.get(url)
        self.assertEqual(response["Content-Type"], "application/json")

    def test_challenge_response_contains_expected_keys(self):
        ''' Test challenge response contains the expecte keys. '''
        drf_client = APIClient(enforce_csrf_checks=True)
        url = reverse('altcha')
        response = drf_client.get(url)
        expected_keys = ["algorithm", "challenge", "maxNumber", "salt", "signature"]
        self.assertEqual(expected_keys, list(response.json().keys()))
