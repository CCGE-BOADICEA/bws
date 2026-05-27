"""
© 2026 University of Cambridge
SPDX-FileCopyrightText: 2026 University of Cambridge
SPDX-License-Identifier: GPL-3.0-or-later
"""

import pytest
from django.test.testcases import TestCase
from django.urls import reverse
from rest_framework import status
from rest_framework.test import APIClient
from unittest.mock import patch, MagicMock


class AltchaWebServices(TestCase):
    ''' Test the altcha webservice '''
                 
    @pytest.mark.req_WS_CORE_200
    def test_challenge_view_returns_200(self):
        ''' Test fetching a new random challenge to be used by the ALTCHA widget. '''
        drf_client = APIClient(enforce_csrf_checks=True)
        url = reverse('altcha')
        response = drf_client.get(url)
        self.assertEqual(response.status_code, status.HTTP_200_OK)

    @pytest.mark.req_WS_CORE_200
    def test_challenge_view_returns_json(self):
        ''' Test challenge response returns JSON. '''
        drf_client = APIClient(enforce_csrf_checks=True)
        url = reverse('altcha')
        response = drf_client.get(url)
        self.assertEqual(response["Content-Type"], "application/json")

    @pytest.mark.req_WS_CORE_200
    def test_challenge_response_contains_expected_keys(self):
        ''' Test challenge response contains the expecte keys. '''
        drf_client = APIClient(enforce_csrf_checks=True)
        url = reverse('altcha')
        response = drf_client.get(url)
        expected_keys = ["algorithm", "challenge", "maxNumber", "salt", "signature"]
        self.assertEqual(expected_keys, list(response.json().keys()))

    @pytest.mark.req_WS_CORE_200
    def test_challenge_response_contains_valid_data_types(self):
        ''' Test challenge response contains fields with valid data types. '''
        drf_client = APIClient(enforce_csrf_checks=True)
        url = reverse('altcha')
        response = drf_client.get(url)
        data = response.json()
        
        # Verify each field has the expected type
        self.assertIsInstance(data.get('algorithm'), str)
        self.assertIsInstance(data.get('challenge'), str)
        self.assertIsInstance(data.get('maxNumber'), int)
        self.assertIsInstance(data.get('salt'), str)
        self.assertIsInstance(data.get('signature'), str)

    @pytest.mark.req_WS_CORE_200
    def test_challenge_response_has_non_empty_values(self):
        ''' Test challenge response contains non-empty values. '''
        drf_client = APIClient(enforce_csrf_checks=True)
        url = reverse('altcha')
        response = drf_client.get(url)
        data = response.json()
        
        # Verify fields are not empty
        self.assertTrue(len(data.get('algorithm', '')) > 0)
        self.assertTrue(len(data.get('challenge', '')) > 0)
        self.assertTrue(data.get('maxNumber', 0) > 0)
        self.assertTrue(len(data.get('salt', '')) > 0)
        self.assertTrue(len(data.get('signature', '')) > 0)

    @pytest.mark.req_WS_CORE_200
    def test_multiple_challenges_have_different_values(self):
        ''' Test that multiple challenge requests return different values. '''
        drf_client = APIClient(enforce_csrf_checks=True)
        url = reverse('altcha')
        
        response1 = drf_client.get(url)
        response2 = drf_client.get(url)
        
        data1 = response1.json()
        data2 = response2.json()
        
        # Challenges should be different (random)
        self.assertNotEqual(data1.get('challenge'), data2.get('challenge'))
        # Salt should be different (random)
        self.assertNotEqual(data1.get('salt'), data2.get('salt'))

    @pytest.mark.req_WS_CORE_200
    def test_challenge_view_only_accepts_get_requests(self):
        ''' Test that challenge endpoint only accepts GET requests '''
        drf_client = APIClient(enforce_csrf_checks=False)
        url = reverse('altcha')
        
        # POST should not be allowed
        response = drf_client.post(url)
        self.assertEqual(response.status_code, status.HTTP_405_METHOD_NOT_ALLOWED)
        
        # PUT should not be allowed
        response = drf_client.put(url)
        self.assertEqual(response.status_code, status.HTTP_405_METHOD_NOT_ALLOWED)
        
        # DELETE should not be allowed
        response = drf_client.delete(url)
        self.assertEqual(response.status_code, status.HTTP_405_METHOD_NOT_ALLOWED)

    @pytest.mark.req_WS_CORE_200
    def test_challenge_view_handles_rapid_requests(self):
        ''' Test that challenge endpoint handles multiple rapid requests. '''
        drf_client = APIClient(enforce_csrf_checks=True)
        url = reverse('altcha')
        
        # Make 10 rapid requests
        responses = []
        for _ in range(10):
            response = drf_client.get(url)
            responses.append(response.status_code)
        
        # All should succeed
        for status_code in responses:
            self.assertEqual(status_code, status.HTTP_200_OK)

    @pytest.mark.req_WS_CORE_200
    def test_challenge_response_algorithm_field_valid(self):
        ''' Test that the algorithm field contains a valid value. '''
        drf_client = APIClient(enforce_csrf_checks=True)
        url = reverse('altcha')
        response = drf_client.get(url)
        data = response.json()
        
        # Algorithm should be SHA-256 or similar
        self.assertIn(data.get('algorithm'), ['SHA-256', 'sha256', 'SHA256'])

    @pytest.mark.req_WS_CORE_200
    @patch('bws.altcha.create_challenge')
    def test_challenge_view_handles_exception(self, mock_create_challenge):
        ''' Test that challenge endpoint handles exceptions gracefully. '''
        mock_create_challenge.side_effect = Exception("Test error")
        
        drf_client = APIClient(enforce_csrf_checks=True)
        url = reverse('altcha')
        response = drf_client.get(url)
        
        # Should return 500 error
        self.assertEqual(response.status_code, status.HTTP_500_INTERNAL_SERVER_ERROR)
        # Should contain error message
        self.assertIn('error', response.json())
