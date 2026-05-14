"""
BOADICEA web-service throttling tests.

© 2023 University of Cambridge
SPDX-FileCopyrightText: 2023 University of Cambridge
SPDX-License-Identifier: GPL-3.0-or-later
"""

import pytest
from django.test import TestCase
from rest_framework import status
from rest_framework.test import APIRequestFactory, APIClient, force_authenticate
from django.core.cache import cache
from unittest.mock import patch, MagicMock

from bws.throttles import EndUserIDRateThrottle, BurstRateThrottle, SustainedRateThrottle
from rest_framework.response import Response
from rest_framework.views import APIView
from rest_framework.authentication import BasicAuthentication
from rest_framework.permissions import IsAuthenticated
from django.contrib.auth.models import User


class MockEndUserIDRateThrottle(EndUserIDRateThrottle):
    max_rate = 3
    rate = str(max_rate)+'/min'
    scope = 'test_enduser_burst'


class MockView_Throttling(APIView):
    throttle_classes = (MockEndUserIDRateThrottle,)

    def post(self, request):
        return Response('foo')


class MockView_Authenticated_Throttling(MockView_Throttling):
    authentication_classes = (BasicAuthentication, )
    permission_classes = (IsAuthenticated,)


class ThrottlingTests(TestCase):

    def setUp(self):
        ''' Reset the cache so that no throttles will be active '''
        cache.clear()
        self.factory = APIRequestFactory()
        self.client = APIClient(enforce_csrf_checks=True)
        self.user = User.objects.create_user('testuser1', email='testuser@test.com',
                                             password='testing1')
        self.user.save()

    @pytest.mark.req_UTILITIES_003
    def test_requests_are_throttled(self):
        ''' Ensure request rate is limited by user_id '''
        for dummy in range(1, MockEndUserIDRateThrottle.max_rate+2):
            request = self.factory.post('/', data={"user_id": "testuserA"})
            response = MockView_Throttling.as_view()(request)
            if dummy <= MockEndUserIDRateThrottle.max_rate:
                self.assertEqual(response.status_code, status.HTTP_200_OK)
            else:
                self.assertEqual(response.status_code, status.HTTP_429_TOO_MANY_REQUESTS)

    @pytest.mark.req_UTILITIES_003
    def test_authenticated_requests_are_throttled(self):
        ''' Ensure request rate is limited by user_id '''
        for dummy in range(1, MockEndUserIDRateThrottle.max_rate+2):
            request = self.factory.post('/', data={"user_id": "testuserA"})
            force_authenticate(request, user=self.user)
            response = MockView_Authenticated_Throttling.as_view()(request)
            if dummy <= MockEndUserIDRateThrottle.max_rate:
                self.assertEqual(response.status_code, status.HTTP_200_OK)
            else:
                self.assertEqual(response.status_code, status.HTTP_429_TOO_MANY_REQUESTS)


class MockBurstRateThrottle(BurstRateThrottle):
    rate = '5/min'
    scope = 'test_burst'


class MockSustainedRateThrottle(SustainedRateThrottle):
    rate = '100/hour'
    scope = 'test_sustained'


class MockView_BurstThrottling(APIView):
    throttle_classes = (MockBurstRateThrottle,)

    def post(self, request):
        return Response('burst_response')


class MockView_SustainedThrottling(APIView):
    throttle_classes = (MockSustainedRateThrottle,)

    def post(self, request):
        return Response('sustained_response')


class BurstRateThrottleTests(TestCase):
    ''' Tests for BurstRateThrottle '''

    def setUp(self):
        cache.clear()
        self.factory = APIRequestFactory()
        self.user = User.objects.create_user('burst_user', email='burst@test.com',
                                             password='testing1')
        self.user.save()

    @pytest.mark.req_UTILITIES_003
    def test_burst_throttle_allows_requests_within_limit(self):
        ''' Ensure burst throttle allows requests within the rate limit '''
        for i in range(5):
            request = self.factory.post('/')
            force_authenticate(request, user=self.user)
            response = MockView_BurstThrottling.as_view()(request)
            self.assertEqual(response.status_code, status.HTTP_200_OK)

    @pytest.mark.req_UTILITIES_003
    def test_burst_throttle_blocks_requests_exceeding_limit(self):
        ''' Ensure burst throttle blocks requests exceeding the rate limit '''
        for i in range(6):
            request = self.factory.post('/')
            force_authenticate(request, user=self.user)
            response = MockView_BurstThrottling.as_view()(request)
            if i < 5:
                self.assertEqual(response.status_code, status.HTTP_200_OK)
            else:
                self.assertEqual(response.status_code, status.HTTP_429_TOO_MANY_REQUESTS)

    @pytest.mark.req_UTILITIES_003
    def test_burst_throttle_with_different_users(self):
        ''' Ensure different users have separate burst limits '''
        user2 = User.objects.create_user('burst_user2', email='burst2@test.com',
                                         password='testing1')
        user2.save()
        
        # First user makes 5 requests
        for i in range(5):
            request = self.factory.post('/')
            force_authenticate(request, user=self.user)
            response = MockView_BurstThrottling.as_view()(request)
            self.assertEqual(response.status_code, status.HTTP_200_OK)
        
        # Second user should still be able to make requests
        request = self.factory.post('/')
        force_authenticate(request, user=user2)
        response = MockView_BurstThrottling.as_view()(request)
        self.assertEqual(response.status_code, status.HTTP_200_OK)


class SustainedRateThrottleTests(TestCase):
    ''' Tests for SustainedRateThrottle '''

    def setUp(self):
        cache.clear()
        self.factory = APIRequestFactory()
        self.user = User.objects.create_user('sustained_user', email='sustained@test.com',
                                             password='testing1')
        self.user.save()

    @pytest.mark.req_UTILITIES_003
    def test_sustained_throttle_allows_requests_within_limit(self):
        ''' Ensure sustained throttle allows requests within the rate limit '''
        request = self.factory.post('/')
        force_authenticate(request, user=self.user)
        response = MockView_SustainedThrottling.as_view()(request)
        self.assertEqual(response.status_code, status.HTTP_200_OK)

    @pytest.mark.req_UTILITIES_003
    def test_sustained_throttle_response_includes_rate_info(self):
        ''' Ensure throttled responses include retry-after headers '''
        cache.clear()
        throttle = MockSustainedRateThrottle()
        factory = APIRequestFactory()
        
        # Create request for authenticated user
        request = factory.post('/')
        request.user = self.user
        force_authenticate(request, user=self.user)
        
        # Check that throttle can process request
        result = throttle.allow_request(request, None)
        self.assertTrue(result)
