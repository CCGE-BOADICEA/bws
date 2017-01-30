from django.test import TestCase
from bws.risk_factors import RiskFactors
from boadicea.exceptions import RiskFactorError
from rest_framework.test import APIRequestFactory, APIClient
from django.contrib.auth.models import User
from boadicea_auth.models import UserDetails
from django.core.urlresolvers import reverse
from rest_framework.authtoken.models import Token
from rest_framework import status
import json
from django.utils.encoding import force_text


class RiskFactorsTests(TestCase):

    def test_round_trip(self):
        ''' Test encoding of risk categories and decoding returns same risk categories. '''
        category1 = list(RiskFactors.categories.values())
        category2 = RiskFactors.decode(RiskFactors.encode(category1))
        self.assertListEqual(category1, category2, "round trip encode/decode")

    def test_wrong_no_risks(self):
        ''' Test that an error is raised when the wrong number of risks is submitted '''
        self.assertRaises(RiskFactorError, RiskFactors.encode, [7, 4, 4, 3, 4])


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
        self.token = Token.objects.create(user=self.user)
        self.token.save()
        self.url = reverse('risk_factors')

    def test_risk_factors(self):
        ''' Test POSTing to the risk factors using token authentication. '''
        data = RiskFactors.categories
        self.client.credentials(HTTP_AUTHORIZATION='Token ' + self.token.key)
        response = self.client.post(self.url, data, format='multipart',
                                    HTTP_ACCEPT="application/json")
        self.assertContains(response, "factor")
        content = json.loads(force_text(response.content))
        factor = RiskFactors.encode(list(data.values()))
        self.assertEqual(content['factor'], factor, "compare result web-service to direct function call")

    def test_risk_factors_bad_request(self):
        ''' Test POSTing to the risk factors with an out of range category number. '''
        ncat = RiskFactors.categories.get('menarche_age')+1
        self.client.credentials(HTTP_AUTHORIZATION='Token ' + self.token.key)
        response = self.client.post(self.url, {'menarche_age': ncat}, format='multipart',
                                    HTTP_ACCEPT="application/json")
        self.assertContains(response, "Ensure this value is less than or equal to ",
                            status_code=status.HTTP_400_BAD_REQUEST)
