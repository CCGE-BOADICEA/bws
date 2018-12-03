""" BOADICEA web-service throttling tests.  """

from django.test import TestCase
from rest_framework import status
from rest_framework.test import APIRequestFactory, APIClient, force_authenticate
from django.core.cache import cache

from bws.throttles import EndUserIDRateThrottle
from rest_framework.response import Response
from rest_framework.views import APIView
from rest_framework.authentication import BasicAuthentication
from rest_framework.permissions import IsAuthenticated
from django.contrib.auth.models import User


class TestEndUserIDRateThrottle(EndUserIDRateThrottle):
    max_rate = 3
    rate = str(max_rate)+'/min'
    scope = 'test_enduser_burst'


class MockView_Throttling(APIView):
    throttle_classes = (TestEndUserIDRateThrottle,)

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

    def test_requests_are_throttled(self):
        ''' Ensure request rate is limited by user_id '''
        for dummy in range(1, TestEndUserIDRateThrottle.max_rate+2):
            request = self.factory.post('/', data={"user_id": "testuserA"})
            response = MockView_Throttling.as_view()(request)
            if dummy <= TestEndUserIDRateThrottle.max_rate:
                self.assertEqual(response.status_code, status.HTTP_200_OK)
            else:
                self.assertEqual(response.status_code, status.HTTP_429_TOO_MANY_REQUESTS)

    def test_authenticated_requests_are_throttled(self):
        ''' Ensure request rate is limited by user_id '''
        for dummy in range(1, TestEndUserIDRateThrottle.max_rate+2):
            request = self.factory.post('/', data={"user_id": "testuserA"})
            force_authenticate(request, user=self.user)
            response = MockView_Authenticated_Throttling.as_view()(request)
            if dummy <= TestEndUserIDRateThrottle.max_rate:
                self.assertEqual(response.status_code, status.HTTP_200_OK)
            else:
                self.assertEqual(response.status_code, status.HTTP_429_TOO_MANY_REQUESTS)
