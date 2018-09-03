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
from bws.risk_factors import bc, oc
from bws.risk_factors.oc import OCRiskFactors


class RiskFactorsCategoryTests(TestCase):

    def test_get_Parity_category(self):
        ''' Given a parity value check the category is correctly assigned. '''
        self.assertEqual(bc.Parity.get_category(3), 4)
        self.assertEqual(oc.Parity.get_category(3), 3)

    def test_get_MHT_category(self):
        ''' Given a MHT value check the category is correctly assigned. '''
        self.assertEqual(bc.MHT.get_category("N"), 1)
        self.assertEqual(oc.MHT.get_category("N"), 1)

        self.assertEqual(bc.MHT.get_category("never/former"), 1)
        self.assertEqual(oc.MHT.get_category("ever"), 2)

    def test_get_BMI_category(self):
        ''' Given a BMI value check the category is correctly assigned. '''
        self.assertEqual(bc.BMI.get_category(" 22 "), 2)
        self.assertEqual(oc.BMI.get_category(22), 1)
        self.assertEqual(bc.BMI.get_category(25), 3)
        self.assertEqual(oc.BMI.get_category(25), 3)
        self.assertEqual(bc.BMI.get_category(30), 4)
        self.assertEqual(oc.BMI.get_category(30), 5)

    def test_get_Height_category(self):
        ''' Given a height value check the category is correctly assigned. '''
        self.assertEqual(bc.Height.get_category(" 152.4 "), 2)
        self.assertEqual(oc.Height.get_category(" 152.4 "), 1)
        self.assertEqual(bc.Height.get_category(173), 4)
        self.assertEqual(oc.Height.get_category(173), 5)

    def test_get_OralContraception_category(self):
        ''' Given a Oral Contraception value check the category is correctly assigned. '''
        self.assertEqual(bc.OralContraception.get_category('-'), 0)
        self.assertEqual(bc.OralContraception.get_category('N'), 1)
        self.assertEqual(bc.OralContraception.get_category('Never'), 1)
        self.assertEqual(bc.OralContraception.get_category('F'), 2)
        self.assertEqual(bc.OralContraception.get_category('C:4'), 3)

        self.assertEqual(oc.OralContraception.get_category('-'), 0)
        self.assertEqual(oc.OralContraception.get_category('N'), 1)
        self.assertEqual(oc.OralContraception.get_category('Never'), 1)
        self.assertEqual(oc.OralContraception.get_category('C:4'), 2)
        self.assertEqual(oc.OralContraception.get_category('C:5'), 3)
        self.assertEqual(oc.OralContraception.get_category('C'), 0)

    def test_get_Endometriosis_category(self):
        ''' Given a endometriosis value check the category is correctly assigned. '''
        self.assertEqual(oc.Endometriosis.get_category('NA'), 0)
        self.assertEqual(oc.Endometriosis.get_category('YES '), 2)
        self.assertEqual(oc.Endometriosis.get_category('n '), 1)

    def test_get_TubalLigation_category(self):
        ''' Given a Tubal Ligation value check the category is correctly assigned. '''
        self.assertEqual(oc.TubalLigation.get_category('-'), 0)
        self.assertEqual(oc.TubalLigation.get_category('no'), 1)
        self.assertEqual(oc.TubalLigation.get_category(' Y'), 2)

    def test_get_FirstLiveBirth_category(self):
        ''' Given a First Live Birth value check the category is correctly assigned. '''
        self.assertEqual(bc.AgeOfFirstLiveBirth.get_category('-'), 0)
        self.assertEqual(bc.AgeOfFirstLiveBirth.get_category(19), 1)
        self.assertEqual(bc.AgeOfFirstLiveBirth.get_category(23), 2)
        self.assertEqual(bc.AgeOfFirstLiveBirth.get_category(32), 4)

    def test_get_Alcohol_category(self):
        ''' Given a Alcohol Intake value check the category is correctly assigned. '''
        self.assertEqual(bc.AlcoholIntake.get_category('NA'), 0)
        self.assertEqual(bc.AlcoholIntake.get_category("0 "), 1)
        self.assertEqual(bc.AlcoholIntake.get_category(3.5), 2)
        self.assertEqual(bc.AlcoholIntake.get_category(33.1), 5)
        self.assertEqual(bc.AlcoholIntake.get_category(45), 7)

    def test_get_AgeOfMenopause_category(self):
        ''' Given an Age Of Menopause value check the category is correctly assigned. '''
        self.assertEqual(bc.AgeOfMenopause.get_category('-'), 0)
        self.assertEqual(bc.AgeOfMenopause.get_category('39'), 1)
        self.assertEqual(bc.AgeOfMenopause.get_category(55), 5)

    def test_get_MammographicDensity_category(self):
        ''' Given a Mammographic Density value check the category is correctly assigned. '''
        self.assertEqual(bc.MammographicDensity.get_category('-'), 0)
        self.assertEqual(bc.MammographicDensity.get_category('3'), 3)
        self.assertEqual(bc.MammographicDensity.get_category("BI-RADS 3"), 3)

    def test_risk_factor_code(self):
        '''
        Test the risk factor code generated
        RFC =    Parity category +
                 Oral Contraception category * 4 +
                 MHT category * 16 +
                 Tubal Ligation category * 48 +
                 Endometriosis category * 144 +
                 BMI category * 432 +
                 Height category * 2592
        '''
        oc_risk_categories = [0 for _k in OCRiskFactors.categories.keys()]

        oc_risk_categories[0] = oc.Parity.get_category('2')
        rfc = 3
        self.assertEqual(OCRiskFactors.encode(oc_risk_categories), rfc)

        oc_risk_categories[1] = oc.OralContraception.get_category('C:4')
        rfc += 2*4
        self.assertEqual(OCRiskFactors.encode(oc_risk_categories), rfc)

        oc_risk_categories[2] = oc.MHT.get_category('C')
        rfc += 2*16
        self.assertEqual(OCRiskFactors.encode(oc_risk_categories), rfc)

        oc_risk_categories[3] = oc.TubalLigation.get_category('no')
        rfc += 1*48
        self.assertEqual(OCRiskFactors.encode(oc_risk_categories), rfc)

        oc_risk_categories[4] = oc.Endometriosis.get_category('yes')
        rfc += 2*144
        self.assertEqual(OCRiskFactors.encode(oc_risk_categories), rfc)

        oc_risk_categories[5] = oc.BMI.get_category(25)
        rfc += 3*432
        self.assertEqual(OCRiskFactors.encode(oc_risk_categories), rfc)

        oc_risk_categories[6] = oc.Height.get_category(153)
        rfc += 2*2592
        self.assertEqual(OCRiskFactors.encode(oc_risk_categories), rfc)


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
