import json

from django.contrib.auth.models import User, Permission
from django.core.urlresolvers import reverse
from django.test import TestCase
from django.utils.encoding import force_text
from rest_framework import status
from rest_framework.authtoken.models import Token
from rest_framework.test import APIRequestFactory, APIClient

from boadicea_auth.models import UserDetails
from bws.exceptions import RiskFactorError
from bws.risk_factors.bc import BCRiskFactors


class RiskFactorsTests(TestCase):

    def test_round_trip(self):
        ''' Test encoding of risk categories and decoding returns same risk categories. '''
        category1 = list(BCRiskFactors.categories.values())
        category2 = BCRiskFactors.decode(BCRiskFactors.encode(category1))
        self.assertListEqual(category1, category2, "round trip encode/decode")

    def test_wrong_no_risks(self):
        ''' Test that an error is raised when the wrong number of risks is submitted '''
        self.assertRaises(RiskFactorError, BCRiskFactors.encode, [7, 4, 4, 3, 4])

    def test_bounds_exceeded_u_encoding(self):
        ''' Test that an error is raised when a risk is above bounds - encoding '''
        category = list(BCRiskFactors.categories.values())
        category[0] += 1
        self.assertRaises(RiskFactorError, BCRiskFactors.encode, category)

    def test_bounds_exceeded_l_encoding(self):
        ''' Test that an error is raised when a risk is below bounds - encoding '''
        category = list(BCRiskFactors.categories.values())
        category[0] = -1
        self.assertRaises(RiskFactorError, BCRiskFactors.encode, category)

    def test_non_numeric_encoding(self):
        ''' Test that an error is raised when passed a non numeric argument - encoding '''
        category = list(BCRiskFactors.categories.values())
        category[0] = 'a'
        self.assertRaises(RiskFactorError, BCRiskFactors.encode, category)

    def test_bounds_exceeded_u_decoding(self):
        ''' Test that an error is raised when a risk is above bounds - decoding '''
        max_plus_1 = BCRiskFactors.encode(list(BCRiskFactors.categories.values())) + 1
        self.assertRaises(RiskFactorError, BCRiskFactors.decode, max_plus_1)

    def test_bounds_exceeded_l_decoding(self):
        ''' Test that an error is raised when a risk is below bounds - decoding '''
        self.assertRaises(RiskFactorError, BCRiskFactors.decode, -1)

    def test_non_numeric_decoding(self):
        ''' Test that an error is raised when passed a non numeric argument - decoding '''
        self.assertRaises(RiskFactorError, BCRiskFactors.decode, 'a')


class RiskFactorsWebServices(TestCase):
    ''' Test the risk factors webservice '''

    def setUp(self):
        ''' Create a user and set up the test client. '''
        self.factory = APIRequestFactory()
        self.client = APIClient(enforce_csrf_checks=True)
        self.user = User.objects.create_user('testuser', email='testuser@test.com',
                                             password='testing')
        # add user details
        UserDetails.objects.create(user=self.user, job_title=UserDetails.CGEN,
                                   country='UK')
        self.user.save()
        self.permission = Permission.objects.get(name='Can risk')
        self.user.user_permissions.add(self.permission)
        self.token = Token.objects.create(user=self.user)
        self.token.save()
        self.url = reverse('internal_ws:risk_factors')

    def test_risk_factors(self):
        ''' Test POSTing to the risk factors using token authentication. '''
        data = BCRiskFactors.categories
        self.client.credentials(HTTP_AUTHORIZATION='Token ' + self.token.key)
        response = self.client.post(self.url, data, format='multipart',
                                    HTTP_ACCEPT="application/json")
        self.assertContains(response, "factor")
        content = json.loads(force_text(response.content))
        factor = BCRiskFactors.encode(list(data.values()))
        self.assertEqual(content['factor'], factor, "compare result web-service to direct function call")

    def test_risk_factors_bad_request(self):
        ''' Test POSTing to the risk factors with an out of range category number. '''
        ncat = BCRiskFactors.categories.get('menarche_age')+1
        self.client.credentials(HTTP_AUTHORIZATION='Token ' + self.token.key)
        response = self.client.post(self.url, {'menarche_age': ncat}, format='multipart',
                                    HTTP_ACCEPT="application/json")
        self.assertContains(response, "Ensure this value is less than or equal to ",
                            status_code=status.HTTP_400_BAD_REQUEST)

    def test_risk_factors_permissions(self):
        ''' Test POSTing to the risk factors without the correct permission. '''
        self.user.user_permissions.remove(self.permission)
        ncat = BCRiskFactors.categories.get('menarche_age')
        self.client.credentials(HTTP_AUTHORIZATION='Token ' + self.token.key)
        response = self.client.post(self.url, {'menarche_age': ncat}, format='multipart',
                                    HTTP_ACCEPT="application/json")
        self.assertContains(response, "Cancer risk factor permission not granted",
                            status_code=status.HTTP_403_FORBIDDEN)
