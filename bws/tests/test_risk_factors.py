"""
© 2022 Cambridge University
SPDX-FileCopyrightText: 2022 Cambridge University
SPDX-License-Identifier: GPL-3.0-or-later
"""
from bws.exceptions import RiskFactorError
from bws.risk_factors import bc, oc
from bws.risk_factors.bc import BCRiskFactors
from bws.risk_factors.oc import OCRiskFactors
from coreapi.compat import force_text
from django.contrib.auth.models import User
from django.test import TestCase
from django.urls import reverse
from rest_framework import status
from rest_framework.authtoken.models import Token
from rest_framework.test import APIClient
import json
import os
from bws.risk_factors.mdensity import Birads, Stratus, Volpara
from bws.risk_factors.ethnicity import UKBioBankEthnicty


class UKBioBankEthnictyTests(TestCase):
    
    def test_white(self):
        ethnicity = UKBioBankEthnicty.factory("White;English/Welsh/Scottish/Northern Irish/British")
        self.assertEqual(ethnicity.get_filename(), "UK-european.nml")

    def test_chinese(self):    
        ethnicity = UKBioBankEthnicty.factory("Asian or Asian British;Chinese")
        self.assertEqual(ethnicity.ethnicity, "chinese")
        
    def test_asian(self):    
        ethnicity = UKBioBankEthnicty.factory("Asian or Asian British;Indian;")
        self.assertEqual(ethnicity.ethnicity, "asian")


class MammographicDensityTests(TestCase):
    '''
    Encoding for the pedigree (Fortran) file:
    If employing BIRAD, the value should be an integer between 1 and 4. If employing
    continuous methods, the value should be a real number in the form N.xxxx. 
    In this case, N refers to the method (10=Stratus, 20= Volpara) whereas xxxxx
    is the mammographic density percentage. 
    Example 1: MD = 42.42% measured with Volpara should be coded as “20.42420”
    Example 2: MD = category 3 of Birads should be coded as “3”
    '''
    
    def test_get_Birads_category(self):
        ''' Given a Birads value check the category is correctly assigned. '''
        self.assertEqual(Birads.get_category('-'), 0)
        self.assertEqual(Birads.get_category('3'), 3)
        self.assertEqual(Birads.get_category("BI-RADS 3"), 3)
        self.assertEqual(Birads.get_category('c'), 3)
        self.assertEqual(Birads.get_category("BI-RADS c"), 3)
        self.assertEqual(Birads.get_category('a'), 1)

    def test_Stratus_pedigree(self):
        ''' Given a Stratus value check the pedigree encoding and display string. '''
        stratus = Stratus("10.67")
        self.assertEqual(stratus.get_pedigree_str(), "10.10670")
        self.assertEqual(stratus.get_display_str(), "Stratus 10.67")
        stratus = Stratus("0.678999")
        self.assertEqual(stratus.get_pedigree_str(), "10.00679")
        self.assertEqual(stratus.get_display_str(), "Stratus 0.678999")
        stratus = Stratus("88")
        self.assertEqual(stratus.get_pedigree_str(), "10.88000")
        self.assertEqual(stratus.get_display_str(), "Stratus 88")

    def test_Volpara_pedigree(self):
        ''' Given a Volpara value check the pedigree encoding and display string. '''
        vol = Volpara("12.67333333")
        self.assertEqual(vol.get_pedigree_str(), "20.12673")
        self.assertEqual(vol.get_display_str(), "Volpara 12.67333333")
        vol = Volpara("1.678999")
        self.assertEqual(vol.get_pedigree_str(), "20.01679")
        self.assertEqual(vol.get_display_str(), "Volpara 1.678999")
        vol = Volpara("5")
        self.assertEqual(vol.get_pedigree_str(), "20.05000")
        self.assertEqual(vol.get_display_str(), "Volpara 5")


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

        # bc_risk_categories[9] = bc.Height.get_category('174.21')
        # rfc += 5*3840000
        self.assertEqual(BCRiskFactors.encode(bc_risk_categories), rfc)

    def test_OC_risk_factor_code(self):
        '''
        Test the ovarian cancer risk factor code generated
        RFC =    Parity category +
                 Oral Contraception category * 4 +
                 MHT category * 24 +
                 Tubal Ligation category * 72 +
                 Endometriosis category * 216 +
                 BMI category * 648
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


class WSRiskFactors(TestCase):
    ''' Test the risk factors webservice '''
    TEST_BASE_DIR = os.path.dirname(os.path.dirname(__file__))
    TEST_DATA_DIR = os.path.join(TEST_BASE_DIR, 'tests', 'data')

    @classmethod
    def setUpClass(cls):
        ''' Create a user, token and url. '''
        super(WSRiskFactors, cls).setUpClass()

        cls.user = User.objects.create_user('testuser', email='testuser@test.com',
                                            password='testing')
        # add user details
        # UserDetails.objects.create(user=cls.user, job_title=UserDetails.CGEN,
        #                           country='UK')
        # cls.user.user_permissions.add(Permission.objects.get(name='Can risk'))
        # cls.user.save()
        cls.token = Token.objects.create(user=cls.user)
        cls.token.save()
        cls.bws_url = reverse('bws')
        cls.ows_url = reverse('ows')

    def setUp(self):
        ''' Set up test client and pedigree data. '''
        self.client = APIClient(enforce_csrf_checks=True)
        with open(os.path.join(WSRiskFactors.TEST_DATA_DIR, "d1.canrisk2"), 'r') as f2:
            self.pedigree_data = f2.read()

    def ws_risk_factor(self, url, bmi, cancer='breast'):
        ''' Test affect of including the risk factors. '''
        data = {'mut_freq': 'UK', 'cancer_rates': 'UK',
                'pedigree_data': self.pedigree_data,
                'user_id': 'test_XXX'}
        self.client.credentials(HTTP_AUTHORIZATION='Token ' + WSRiskFactors.token.key)
        response = self.client.post(url, data, format='multipart',
                                    HTTP_ACCEPT="application/json")
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        cancer_risks1 = json.loads(force_text(response.content))['pedigree_result'][0]['cancer_risks']

        # add BMI
        data['pedigree_data'] = self.pedigree_data.replace("##CanRisk 2.0", f"##CanRisk 2.0\n##BMI={bmi}")
        response = self.client.post(url, data, format='multipart',
                                    HTTP_ACCEPT="application/json")
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        cancer_risks2 = json.loads(force_text(response.content))['pedigree_result'][0]['cancer_risks']
        self.assertLess(cancer_risks1[0][f'{cancer} cancer risk']['decimal'],
                        cancer_risks2[0][f'{cancer} cancer risk']['decimal'])

    def test_bws_risk_factor(self):
        ''' Test affect of including the risk factors in BWS. '''
        self.ws_risk_factor(WSRiskFactors.bws_url, bc.BMI.get_category(30)*3200)

    def test_ows_risk_factor(self):
        ''' Test affect of including the risk factors in oWS. '''
        self.ws_risk_factor(WSRiskFactors.ows_url, oc.BMI.get_category(30)*3200, cancer='ovarian')

    def ws_prs(self, url, cancer='breast'):
        ''' Test affect of including the PRS. '''
        data = {'mut_freq': 'UK', 'cancer_rates': 'UK',
                'pedigree_data': self.pedigree_data,
                'user_id': 'test_XXX'}
        self.client.credentials(HTTP_AUTHORIZATION='Token ' + WSRiskFactors.token.key)
        response = self.client.post(url, data, format='multipart', HTTP_ACCEPT="application/json")
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        cancer_risks1 = json.loads(force_text(response.content))['pedigree_result'][0]['cancer_risks']

        # add PRS
        data['prs'] = json.dumps({'alpha': 0.45, 'zscore': 2.652})
        response = self.client.post(url, data, format='multipart', HTTP_ACCEPT="application/json")
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        cancer_risks2 = json.loads(force_text(response.content))['pedigree_result'][0]['cancer_risks']
        self.assertLess(cancer_risks1[0][f'{cancer} cancer risk']['decimal'],
                        cancer_risks2[0][f'{cancer} cancer risk']['decimal'])

    def test_bws_prs(self):
        ''' Test affect of including the PRS in BWS. '''
        self.ws_prs(WSRiskFactors.bws_url)

    def test_ows_prs(self):
        ''' Test affect of including the PRS in OWS. '''
        self.ws_prs(WSRiskFactors.ows_url, cancer='ovarian')

#    def test_risk_factors_inconsistent(self):
#        ''' Test inconsistent risk factors, e.g. age of first birth specified with parity unobserved. '''
#        data = {'mut_freq': 'UK', 'cancer_rates': 'UK',
#                'pedigree_data': self.pedigree_data,
#                'user_id': 'test_XXX', 'risk_factor_code': BCRiskFactors.encode([0, 0, 1, 0, 0, 0, 0, 0, 0])}
#        self.client.credentials(HTTP_AUTHORIZATION='Token ' + BwsRiskFactors.token.key)
#        BwsRiskFactors.user.user_permissions.add(Permission.objects.get(name='Can risk'))
#        response = self.client.post(BwsRiskFactors.url, data, format='multipart',
#                                    HTTP_ACCEPT="application/json")
#        self.assertEqual(response.status_code, status.HTTP_400_BAD_REQUEST)
#        content = json.loads(force_text(response.content))
#        self.assertTrue('Model Error' in content)
