"""
© 2026 University of Cambridge
SPDX-FileCopyrightText: 2026 University of Cambridge
SPDX-License-Identifier: GPL-3.0-or-later
"""

import time

import pytest
from altcha import Challenge, Payload, solve_challenge, verify_solution
from django.conf import settings
from django.test import override_settings
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
        ''' Test challenge response contains the expected proof-of-work v2 keys. '''
        drf_client = APIClient(enforce_csrf_checks=True)
        url = reverse('altcha')
        response = drf_client.get(url)
        data = response.json()

        self.assertEqual(["parameters", "signature"], list(data.keys()))
        expected_keys = ["algorithm", "cost", "keyLength", "keyPrefix", "nonce", "salt", "expiresAt"]
        self.assertEqual(sorted(expected_keys), sorted(data["parameters"].keys()))

    @pytest.mark.req_WS_CORE_200
    def test_challenge_response_contains_valid_data_types(self):
        ''' Test challenge response contains fields with valid data types. '''
        drf_client = APIClient(enforce_csrf_checks=True)
        url = reverse('altcha')
        response = drf_client.get(url)
        data = response.json()
        params = data["parameters"]

        # Verify each field has the expected type
        self.assertIsInstance(params.get('algorithm'), str)
        self.assertIsInstance(params.get('cost'), int)
        self.assertIsInstance(params.get('keyLength'), int)
        self.assertIsInstance(params.get('keyPrefix'), str)
        self.assertIsInstance(params.get('nonce'), str)
        self.assertIsInstance(params.get('salt'), str)
        self.assertIsInstance(params.get('expiresAt'), int)
        self.assertIsInstance(data.get('signature'), str)

    @pytest.mark.req_WS_CORE_200
    def test_challenge_response_has_non_empty_values(self):
        ''' Test challenge response contains non-empty values. '''
        drf_client = APIClient(enforce_csrf_checks=True)
        url = reverse('altcha')
        response = drf_client.get(url)
        data = response.json()
        params = data["parameters"]

        # Verify fields are not empty
        self.assertTrue(len(params.get('algorithm', '')) > 0)
        self.assertTrue(params.get('cost', 0) > 0)
        self.assertTrue(len(params.get('nonce', '')) > 0)
        self.assertTrue(len(params.get('salt', '')) > 0)
        self.assertTrue(len(data.get('signature', '')) > 0)

    @pytest.mark.req_WS_CORE_200
    def test_challenge_expires_in_the_future(self):
        ''' Test the challenge carries an expiry, so a solution can not be replayed forever. '''
        drf_client = APIClient(enforce_csrf_checks=True)
        url = reverse('altcha')
        response = drf_client.get(url)
        self.assertGreater(response.json()["parameters"]["expiresAt"], int(time.time()))

    @pytest.mark.req_WS_CORE_200
    def test_multiple_challenges_have_different_values(self):
        ''' Test that multiple challenge requests return different values. '''
        drf_client = APIClient(enforce_csrf_checks=True)
        url = reverse('altcha')

        response1 = drf_client.get(url)
        response2 = drf_client.get(url)

        data1 = response1.json()["parameters"]
        data2 = response2.json()["parameters"]

        # Nonces should be different (random)
        self.assertNotEqual(data1.get('nonce'), data2.get('nonce'))
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
        ''' Test that the algorithm field contains the configured key derivation function. '''
        drf_client = APIClient(enforce_csrf_checks=True)
        url = reverse('altcha')
        response = drf_client.get(url)
        data = response.json()

        self.assertEqual(data["parameters"].get('algorithm'), settings.ALTCHA_ALGORITHM)

    @pytest.mark.req_WS_CORE_200
    @override_settings(ALTCHA_COST=10)      # keep the proof-of-work cheap enough to solve here
    def test_challenge_can_be_solved_and_verified(self):
        ''' Test a challenge issued by the endpoint can be solved and then verified, i.e. that
        the endpoint and the form verification (see AltchaFormMixin) agree. '''
        drf_client = APIClient(enforce_csrf_checks=True)
        url = reverse('altcha')
        challenge = Challenge.from_dict(drf_client.get(url).json())

        solution = solve_challenge(challenge)
        self.assertIsNotNone(solution, "challenge could not be solved")

        payload = Payload(challenge=challenge, solution=solution).to_base64()
        result = verify_solution(payload, settings.ALTCHA_HMAC_KEY)
        self.assertTrue(result.verified)
        self.assertFalse(result.expired)

    @pytest.mark.req_WS_CORE_200
    @override_settings(ALTCHA_COST=10)
    def test_solution_rejected_with_wrong_hmac_key(self):
        ''' Test a solution is rejected when verified with a different HMAC key, i.e. a
        challenge that was not issued by this site. '''
        drf_client = APIClient(enforce_csrf_checks=True)
        challenge = Challenge.from_dict(drf_client.get(reverse('altcha')).json())
        payload = Payload(challenge=challenge, solution=solve_challenge(challenge)).to_base64()

        result = verify_solution(payload, 'not-the-hmac-key')
        self.assertFalse(result.verified)
        self.assertTrue(result.invalid_signature)

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
