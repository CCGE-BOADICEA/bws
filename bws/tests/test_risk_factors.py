import json
import os

from django.contrib.auth.models import User, Permission
from django.core.urlresolvers import reverse
from django.test import TestCase
from django.utils.encoding import force_text
from rest_framework import status
from rest_framework.authtoken.models import Token
from rest_framework.test import APIClient

from boadicea_auth.models import UserDetails
from bws.exceptions import RiskFactorError
from bws.risk_factors import bc, oc
from bws.risk_factors.bc import BCRiskFactors
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

        self.assertEqual(oc.BMI.get_category(25), 2)
        self.assertEqual(bc.BMI.get_category(30), 4)
        self.assertEqual(oc.BMI.get_category(30), 3)

    def test_get_Height_category(self):
        ''' Given a height value check the category is correctly assigned. '''
        self.assertEqual(bc.Height.get_category(" 152.91 "), 2)
        self.assertEqual(oc.Height.get_category(" 152.4 "), 1)
        self.assertEqual(bc.Height.get_category(159.65), 3)
        self.assertEqual(bc.Height.get_category(172.68), 4)
        self.assertEqual(oc.Height.get_category(173), 5)

    def test_get_OralContraception_category(self):
        ''' Given a Oral Contraception value check the category is correctly assigned. '''
        self.assertEqual(bc.OralContraception.get_category('-'), 0)
        self.assertEqual(bc.OralContraception.get_category('N'), 1)
        self.assertEqual(bc.OralContraception.get_category('Never'), 1)
        self.assertEqual(bc.OralContraception.get_category('F'), 2)
        self.assertEqual(bc.OralContraception.get_category('C:4'), 3)
        self.assertEqual(bc.OralContraception.get_category('C'), 3)

        self.assertEqual(oc.OralContraception.get_category('-'), 0)
        self.assertEqual(oc.OralContraception.get_category('N'), 1)
        self.assertEqual(oc.OralContraception.get_category('Never'), 1)
        self.assertEqual(oc.OralContraception.get_category('C:4'), 2)
        self.assertEqual(oc.OralContraception.get_category('C:5'), 3)
        self.assertEqual(oc.OralContraception.get_category('C:0.5'), 1)
        self.assertEqual(oc.OralContraception.get_category('C:<1'), 1)
        self.assertEqual(oc.OralContraception.get_category('C'), 0)
        self.assertEqual(oc.OralContraception.get_category('F'), 0)

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
        self.assertEqual(bc.AlcoholIntake.get_category(44.4), 6)
        self.assertEqual(bc.AlcoholIntake.get_category(44.5), 6)
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


class RiskFactorsCodeTests(TestCase):

    def test_BC_risk_factor_code(self):
        '''
        Test the breast cancer risk factor code generated
        RFC =   menarch category +
                parity category * 8 +
                first birth category * 40 +
                OC category * 200 +
                MHT category * 800 +
                BMI category * 3200 +
                alcohol category * 16000 +
                menopause category * 128000 +
                density category * 768000 +
                height category * 3840000
        '''
        bc_risk_categories = [0 for _k in BCRiskFactors.categories.keys()]
        bc_risk_categories[0] = bc.MenarcheAge.get_category('12')
        rfc = 3
        self.assertEqual(BCRiskFactors.encode(bc_risk_categories), rfc)

        bc_risk_categories[1] = bc.Parity.get_category('1')
        rfc += 2*8
        self.assertEqual(BCRiskFactors.encode(bc_risk_categories), rfc)

        bc_risk_categories[2] = bc.AgeOfFirstLiveBirth.get_category('31')
        rfc += 4*40
        self.assertEqual(BCRiskFactors.encode(bc_risk_categories), rfc)

        bc_risk_categories[3] = bc.OralContraception.get_category('F:3')
        rfc += 2*200
        self.assertEqual(BCRiskFactors.encode(bc_risk_categories), rfc)

        bc_risk_categories[4] = bc.MHT.get_category('C')
        rfc += 3*800
        self.assertEqual(BCRiskFactors.encode(bc_risk_categories), rfc)

        bc_risk_categories[5] = bc.BMI.get_category(22.3)
        rfc += 2*3200
        self.assertEqual(BCRiskFactors.encode(bc_risk_categories), rfc)

        bc_risk_categories[6] = bc.AlcoholIntake.get_category(0)
        rfc += 1*16000
        self.assertEqual(BCRiskFactors.encode(bc_risk_categories), rfc)

        bc_risk_categories[7] = bc.AgeOfMenopause.get_category(52)
        rfc += 4*128000
        self.assertEqual(BCRiskFactors.encode(bc_risk_categories), rfc)

        bc_risk_categories[8] = bc.MammographicDensity.get_category('1')
        rfc += 1*768000
        self.assertEqual(BCRiskFactors.encode(bc_risk_categories), rfc)

        bc_risk_categories[9] = bc.Height.get_category('174.21')
        rfc += 5*3840000
        self.assertEqual(BCRiskFactors.encode(bc_risk_categories), rfc)

    def test_OC_risk_factor_code(self):
        '''
        Test the ovarian cancer risk factor code generated
        RFC =    Parity category +
                 Oral Contraception category * 4 +
                 MHT category * 24 +
                 Tubal Ligation category * 72 +
                 Endometriosis category * 216 +
                 BMI category * 648 +
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
        rfc += 2*24
        self.assertEqual(OCRiskFactors.encode(oc_risk_categories), rfc)

        oc_risk_categories[3] = oc.TubalLigation.get_category('no')
        rfc += 1*72
        self.assertEqual(OCRiskFactors.encode(oc_risk_categories), rfc)

        oc_risk_categories[4] = oc.Endometriosis.get_category('yes')
        rfc += 2*216
        self.assertEqual(OCRiskFactors.encode(oc_risk_categories), rfc)

        oc_risk_categories[5] = oc.BMI.get_category(25)
        rfc += 2*648
        self.assertEqual(OCRiskFactors.encode(oc_risk_categories), rfc)

        oc_risk_categories[6] = oc.Height.get_category(153)
        rfc += 2*2592
        self.assertEqual(OCRiskFactors.encode(oc_risk_categories), rfc)

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


class BwsRiskFactors(TestCase):
    ''' Test the risk factors webservice '''
    TEST_BASE_DIR = os.path.dirname(os.path.dirname(__file__))
    TEST_DATA_DIR = os.path.join(TEST_BASE_DIR, 'tests', 'data')

    @classmethod
    def setUpClass(cls):
        ''' Create a user, token and url. '''
        super(BwsRiskFactors, cls).setUpClass()

        cls.user = User.objects.create_user('testuser', email='testuser@test.com',
                                            password='testing')
        # add user details
        UserDetails.objects.create(user=cls.user, job_title=UserDetails.CGEN,
                                   country='UK')
        cls.user.save()
        cls.token = Token.objects.create(user=cls.user)
        cls.token.save()
        cls.url = reverse('bws')

    def setUp(self):
        ''' Set up test client and pedigree data. '''
        self.client = APIClient(enforce_csrf_checks=True)
        self.pedigree_data = open(os.path.join(BwsRiskFactors.TEST_DATA_DIR, "pedigree_data.txt"), "r")

    def test_bws_risk_factor(self):
        ''' Test affect of including the risk factors. '''
        data = {'mut_freq': 'UK', 'cancer_rates': 'UK',
                'pedigree_data': self.pedigree_data,
                'user_id': 'test_XXX', 'risk_factor_code': 7}
        self.client.credentials(HTTP_AUTHORIZATION='Token ' + BwsRiskFactors.token.key)
        # no permissions to use the risk factors and so ignored
        response = self.client.post(BwsRiskFactors.url, data, format='multipart',
                                    HTTP_ACCEPT="application/json")
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        cancer_risks1 = json.loads(force_text(response.content))['pedigree_result'][0]['cancer_risks']

        # add permissions to enable use of the risk factors
        data['pedigree_data'] = open(os.path.join(BwsRiskFactors.TEST_DATA_DIR, "pedigree_data.txt"), "r")
        BwsRiskFactors.user.user_permissions.add(Permission.objects.get(name='Can risk'))
        response = self.client.post(BwsRiskFactors.url, data, format='multipart',
                                    HTTP_ACCEPT="application/json")
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        cancer_risks2 = json.loads(force_text(response.content))['pedigree_result'][0]['cancer_risks']
        self.assertLess(cancer_risks2[0]['breast cancer risk']['decimal'],
                        cancer_risks1[0]['breast cancer risk']['decimal'])

    def test_risk_factors_inconsistent(self):
        ''' Test inconsistent risk factors, e.g. age of first birth specified with parity unobserved. '''
        data = {'mut_freq': 'UK', 'cancer_rates': 'UK',
                'pedigree_data': self.pedigree_data,
                'user_id': 'test_XXX', 'risk_factor_code': BCRiskFactors.encode([0, 0, 1, 0, 0, 0, 0, 0, 0, 0])}
        self.client.credentials(HTTP_AUTHORIZATION='Token ' + BwsRiskFactors.token.key)
        BwsRiskFactors.user.user_permissions.add(Permission.objects.get(name='Can risk'))
        response = self.client.post(BwsRiskFactors.url, data, format='multipart',
                                    HTTP_ACCEPT="application/json")
        self.assertEqual(response.status_code, status.HTTP_400_BAD_REQUEST)
        content = json.loads(force_text(response.content))
        self.assertTrue('Model Error' in content)
