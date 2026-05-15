"""
Tests for rest_api.py module.

© 2026 University of Cambridge
SPDX-FileCopyrightText: 2026 University of Cambridge
SPDX-License-Identifier: GPL-3.0-or-later
"""
import pytest
import json
import os
from unittest.mock import MagicMock, patch, mock_open
from django.contrib.auth.models import User, Permission
from django.core.exceptions import PermissionDenied
from django.test import TestCase, RequestFactory
from django.urls import reverse
from django.utils.encoding import force_str
from rest_framework import status
from rest_framework.authtoken.models import Token
from rest_framework.test import APIClient, APIRequestFactory
from rest_framework.exceptions import ValidationError
from types import SimpleNamespace

from bws.rest_api import RequiredAnyPermission, ModelWebServiceMixin, BwsView, OwsView, PwsView, CombineModelResultsView
from bws.exceptions import ModelError, PedigreeError
from bws.serializers import CombinedInputSerializer
from django.conf import settings


class TestRequiredAnyPermission(TestCase):
    """Tests for the RequiredAnyPermission class."""

    def setUp(self):
        self.factory = RequestFactory()
        self.permission = RequiredAnyPermission()

    @pytest.mark.req_WS_CORE_108
    def test_has_permission_with_valid_perm(self):
        """Test that permission is granted when user has one of the required permissions."""
        user = User.objects.create_user('testuser', email='test@example.com', password='test')
        user.user_permissions.add(Permission.objects.get(name='Can risk'))
        user.save()

        request = self.factory.get('/')
        request.user = user

        # Mock view with any_perms
        view = MagicMock()
        view.any_perms = ['boadicea_auth.can_risk', 'boadicea_auth.commercial_api_breast']

        self.assertTrue(self.permission.has_permission(request, view))

    @pytest.mark.req_WS_CORE_108
    def test_has_permission_without_valid_perm(self):
        """Test that permission is denied when user has none of the required permissions."""
        user = User.objects.create_user('testuser2', email='test2@example.com', password='test')
        user.save()

        request = self.factory.get('/')
        request.user = user

        # Mock view with any_perms
        view = MagicMock()
        view.any_perms = ['boadicea_auth.can_risk', 'boadicea_auth.commercial_api_breast']

        with self.assertRaises(PermissionDenied):
            self.permission.has_permission(request, view)

    @pytest.mark.req_WS_CORE_108
    def test_has_permission_with_multiple_valid_perms(self):
        """Test that permission is granted when user has multiple required permissions."""
        user = User.objects.create_user('testuser3', email='test3@example.com', password='test')
        user.user_permissions.add(Permission.objects.get(name='Can risk'))
        user.user_permissions.add(Permission.objects.get(name='Commercial BC webservices'))
        user.save()

        request = self.factory.get('/')
        request.user = user

        # Mock view with any_perms
        view = MagicMock()
        view.any_perms = ['boadicea_auth.can_risk', 'boadicea_auth.commercial_api_breast']

        self.assertTrue(self.permission.has_permission(request, view))


class TestModelWebServiceMixin(TestCase):
    """Tests for the ModelWebServiceMixin class."""

    def setUp(self):
        self.mixin = ModelWebServiceMixin()
        self.factory = RequestFactory()

    @pytest.mark.req_WS_CORE_109
    def test_get_risk_factors_bc_model(self):
        """Test get_risk_factors method for BC model."""
        from bws.risk_factors.bc import BCRiskFactors
        risk_factor_code = 0  # All categories set to 0 (first category)

        result = self.mixin.get_risk_factors(settings.BC_MODEL, risk_factor_code)

        # Should return a dictionary with risk factor names as keys
        self.assertIsInstance(result, dict)
        # Check that we get the expected risk factors for BC model
        expected_factors = ['menarche_age', 'parity', 'age_of_first_live_birth',
                            'oral_contraception', 'mht', 'bmi', 'alcohol_intake',
                            'age_of_menopause']
        for factor in expected_factors:
            self.assertIn(factor, result)

    @pytest.mark.req_WS_CORE_109
    def test_get_risk_factors_oc_model(self):
        """Test get_risk_factors method for OC model."""
        risk_factor_code = 0  # All categories set to 0 (first category)

        result = self.mixin.get_risk_factors(settings.OC_MODEL, risk_factor_code)

        # Should return a dictionary with risk factor names as keys
        self.assertIsInstance(result, dict)
        # Check that we get the expected risk factors for OC model
        expected_factors = ['parity', 'oral_contraception', 'mht',
                            'tubal_ligation', 'endometriosis', 'bmi']
        for factor in expected_factors:
            self.assertIn(factor, result)

    @pytest.mark.req_WS_CORE_109
    def test_get_risk_factors_pc_model(self):
        """Test get_risk_factors method for PC model."""
        risk_factor_code = 0  # All categories set to 0 (first category)

        result = self.mixin.get_risk_factors(settings.PC_MODEL, risk_factor_code)

        # Should return a dictionary with risk factor names as keys
        self.assertIsInstance(result, dict)
        # PC model should have no risk factors
        self.assertEqual(len(result), 0)

    @pytest.mark.req_WS_CORE_109
    def test_get_risk_factors_invalid_model(self):
        """Test get_risk_factors method with invalid model name."""
        invalid_model = {'NAME': 'INVALID'}

        with self.assertRaises(ModelError):
            self.mixin.get_risk_factors(invalid_model, 0)

    @pytest.mark.req_WS_CORE_110
    def test_add_attr_success(self):
        """Test add_attr method when attribute exists."""
        calcs = MagicMock()
        calcs.test_attr = "test_value"

        this_pedigree = {}
        output = {}

        self.mixin.add_attr("test_attr", this_pedigree, calcs, output)
        self.assertEqual(this_pedigree["test_attr"], "test_value")
        self.assertNotIn("warnings", output)

    @pytest.mark.req_WS_CORE_110
    def test_add_attr_missing_creates_warning(self):
        """Test add_attr method when attribute is missing."""
        calcs = SimpleNamespace()

        this_pedigree = {}
        output = {}

        self.mixin.add_attr("missing_attr", this_pedigree, calcs, output)
        self.assertNotIn("missing_attr", this_pedigree)
        self.assertIn("warnings", output)
        self.assertIn("missing_attr not provided", output["warnings"])

    @pytest.mark.req_WS_CORE_110
    def test_add_attr_missing_extends_existing_warnings(self):
        """Test add_attr method when warnings already exist."""
        calcs = SimpleNamespace()

        this_pedigree = {}
        output = {"warnings": ["existing warning"]}

        self.mixin.add_attr("missing_attr", this_pedigree, calcs, output)
        self.assertNotIn("missing_attr", this_pedigree)
        self.assertEqual(len(output["warnings"]), 2)
        self.assertIn("existing warning", output["warnings"])
        self.assertIn("missing_attr not provided", output["warnings"])


class TestCombineModelResultsView(TestCase):
    """Tests for the CombineModelResultsView class."""

    def setUp(self):
        self.factory = RequestFactory()
        self.user = User.objects.create_user('testuser', email='test@example.com', password='test')
        self.user.user_permissions.add(Permission.objects.get(name='Can risk'))
        self.user.save()
        self.token = Token.objects.create(user=self.user)
        self.token.save()
        self.client = APIClient(enforce_csrf_checks=True)
        self.client.credentials(HTTP_AUTHORIZATION='Token ' + self.token.key)
        self.url = reverse('combine')

    @pytest.mark.req_WS_CORE_111
    def test_post_valid_data(self):
        """Test POST method with valid combined results data."""
        # Mock valid data for combining results

        valid_data = {
            "ows_result": {
                "version":"ovarian model 3.2.2, version 0.6.0",
                "timestamp":"2026-05-15T13:54:27.486660+01:00",
                "mutation_frequency":{"UK":{"BRCA1":0.0007947,"BRCA2":0.002576,"RAD51D":0.00035,"RAD51C":0.00035,"BRIP1":0.00071,"PALB2":0.00064}},
                "mutation_sensitivity":{"BRCA1":0.89,"BRCA2":0.96,"RAD51D":0.86,"RAD51C":0.78,"BRIP1":0.95,"PALB2":0.92},
                "cancer_incidence_rates":"UK",
                "warnings":["ten_yr_cancer_risk not provided","baseline_ten_yr_cancer_risk not provided","year of birth and age at last follow up must be specified in order for m21, f21 to be included in a calculation"],
                "pedigree_result":[
                    {"family_id":"XXXX","proband_id":"ch1",
                     "mutation_frequency":{"UK":{"BRCA1":0.0007947,"BRCA2":0.002576,"RAD51D":0.00035,"RAD51C":0.00035,"BRIP1":0.00071,"PALB2":0.00064}},
                     "risk_factors":{"parity":"-","oral_contraception":"-","mht":"-","tubal_ligation":"-","endometriosis":"-","bmi":"-","Height (cm)":"-"},
                     "cancer_risks":[{"age":65,"ovarian cancer risk":{"decimal":0.0004829,"percent":0}},{"age":66,"ovarian cancer risk":{"decimal":0.0009869,"percent":0.1}},{"age":67,"ovarian cancer risk":{"decimal":0.0015119,"percent":0.2}},{"age":68,"ovarian cancer risk":{"decimal":0.0020575,"percent":0.2}},{"age":69,"ovarian cancer risk":{"decimal":0.0026235,"percent":0.3}},{"age":70,"ovarian cancer risk":{"decimal":0.0032095,"percent":0.3}},{"age":74,"ovarian cancer risk":{"decimal":0.0057386,"percent":0.6}},{"age":75,"ovarian cancer risk":{"decimal":0.0064123,"percent":0.6}},{"age":80,"ovarian cancer risk":{"decimal":0.0098928,"percent":1}}],
                     "baseline_cancer_risks":[{"age":65,"ovarian cancer risk":{"decimal":0.0004829,"percent":0}},{"age":66,"ovarian cancer risk":{"decimal":0.0009869,"percent":0.1}},{"age":67,"ovarian cancer risk":{"decimal":0.0015119,"percent":0.2}},{"age":68,"ovarian cancer risk":{"decimal":0.0020575,"percent":0.2}},{"age":69,"ovarian cancer risk":{"decimal":0.0026235,"percent":0.3}},{"age":70,"ovarian cancer risk":{"decimal":0.0032095,"percent":0.3}},{"age":74,"ovarian cancer risk":{"decimal":0.0057386,"percent":0.6}},{"age":75,"ovarian cancer risk":{"decimal":0.0064123,"percent":0.6}},{"age":80,"ovarian cancer risk":{"decimal":0.0098928,"percent":1}}],
                     "lifetime_cancer_risk":[{"age":80,"ovarian cancer risk":{"decimal":0.0176191,"percent":1.8}}],
                     "baseline_lifetime_cancer_risk":[{"age":80,"ovarian cancer risk":{"decimal":0.0176191,"percent":1.8}}],
                     "mutation_probabilties":[{"no mutation":{"decimal":0.9933273,"percent":99.33}},{"BRCA1":{"decimal":0.0005071,"percent":0.05}},{"BRCA2":{"decimal":0.0026639,"percent":0.27}},{"RAD51D":{"decimal":0.0006242,"percent":0.06}},{"RAD51C":{"decimal":0.0006159,"percent":0.06}},{"BRIP1":{"decimal":0.0013898,"percent":0.14}},{"PALB2":{"decimal":0.0008718,"percent":0.09}}]}]},
            "bws_result": {
                "version":"boadicea model 7.3.2, version 0.6.0",
                "timestamp":"2026-05-15T13:54:27.486498+01:00",
                "mutation_frequency":{"UK":{"BRCA1":0.0006394,"BRCA2":0.00102,"PALB2":0.00064,"ATM":0.0018,"CHEK2":0.00373,"BARD1":0.00043,"RAD51C":0.00035,"RAD51D":0.00035}},
                "mutation_sensitivity":{"BRCA1":0.89,"BRCA2":0.96,"PALB2":0.92,"ATM":0.94,"CHEK2":0.98,"BARD1":0.89,"RAD51C":0.78,"RAD51D":0.86},
                "cancer_incidence_rates":"UK",
                "warnings":["year of birth and age at last follow up must be specified in order for m21, f21 to be included in a calculation"],
                "pedigree_result":[
                    {"family_id":"XXXX","proband_id":"ch1",
                     "mutation_frequency":{"UK":{"BRCA1":0.0006394,"BRCA2":0.00102,"PALB2":0.00064,"ATM":0.0018,"CHEK2":0.00373,"BARD1":0.00043,"RAD51C":0.00035,"RAD51D":0.00035}},
                     "risk_factors":{"menarche_age":"-","parity":"-","age_of_first_live_birth":"-","oral_contraception":"-","mht":"-","bmi":"-","alcohol_intake":"-","age_of_menopause":"-","Mammographic Density":"-","Height (cm)":"-"},
                     "cancer_risks":[{"age":65,"breast cancer risk":{"decimal":0.0036971,"percent":0.4}},{"age":66,"breast cancer risk":{"decimal":0.007466,"percent":0.7}},{"age":67,"breast cancer risk":{"decimal":0.0112766,"percent":1.1}},{"age":68,"breast cancer risk":{"decimal":0.0151054,"percent":1.5}},{"age":69,"breast cancer risk":{"decimal":0.0189373,"percent":1.9}},{"age":70,"breast cancer risk":{"decimal":0.0227619,"percent":2.3}},{"age":74,"breast cancer risk":{"decimal":0.0378265,"percent":3.8}},{"age":75,"breast cancer risk":{"decimal":0.0415487,"percent":4.2}},{"age":80,"breast cancer risk":{"decimal":0.0606968,"percent":6.1}}],
                     "baseline_cancer_risks":[{"age":65,"breast cancer risk":{"decimal":0.0036971,"percent":0.4}},{"age":66,"breast cancer risk":{"decimal":0.007466,"percent":0.7}},{"age":67,"breast cancer risk":{"decimal":0.0112766,"percent":1.1}},{"age":68,"breast cancer risk":{"decimal":0.0151054,"percent":1.5}},{"age":69,"breast cancer risk":{"decimal":0.0189373,"percent":1.9}},{"age":70,"breast cancer risk":{"decimal":0.0227619,"percent":2.3}},{"age":74,"breast cancer risk":{"decimal":0.0378265,"percent":3.8}},{"age":75,"breast cancer risk":{"decimal":0.0415487,"percent":4.2}},{"age":80,"breast cancer risk":{"decimal":0.0606968,"percent":6.1}}],
                     "lifetime_cancer_risk":[{"age":80,"breast cancer risk":{"decimal":0.11897,"percent":11.9}}],
                     "baseline_lifetime_cancer_risk":[{"age":80,"breast cancer risk":{"decimal":0.11897,"percent":11.9}}],
                     "ten_yr_cancer_risk":[{"age":50,"breast cancer risk":{"decimal":0.0165779,"percent":1.7}}],
                     "baseline_ten_yr_cancer_risk":[{"age":50,"breast cancer risk":{"decimal":0.0165779,"percent":1.7}}],
                     "mutation_probabilties":[{"no mutation":{"decimal":0.9854749,"percent":98.55}},{"BRCA1":{"decimal":0.000408,"percent":0.04}},{"BRCA2":{"decimal":0.0010558,"percent":0.11}},{"PALB2":{"decimal":0.0008746,"percent":0.09}},{"CHEK2":{"decimal":0.0068254,"percent":0.68}},{"ATM":{"decimal":0.0033317,"percent":0.33}},{"BARD1":{"decimal":0.0007955,"percent":0.08}},{"RAD51C":{"decimal":0.0006131,"percent":0.06}},{"RAD51D":{"decimal":0.000621,"percent":0.06}}]}]}
        }

        response = self.client.post(self.url, valid_data, format='json')
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        # Should return HTML template response
        self.assertEqual(response.accepted_renderer.format, 'html')

    @pytest.mark.req_WS_CORE_111
    def test_post_invalid_data(self):
        """Test POST method with invalid data."""
        invalid_data = {"invalid": "data"}
        response = self.client.post(self.url, invalid_data, format='json')
        self.assertEqual(response.status_code, status.HTTP_400_BAD_REQUEST)

    @pytest.mark.req_WS_CORE_111
    def test_post_unauthenticated(self):
        """Test POST method without authentication."""
        unauth_client = APIClient()
        valid_data = {
            "bws_result": {},
            "ows_result": {}
        }

        response = unauth_client.post(self.url, valid_data, format='json')
        self.assertEqual(response.status_code, status.HTTP_403_FORBIDDEN)


class TestModelWebServiceMixinIntegration(TestCase):
    """Integration tests for ModelWebServiceMixin error handling."""

    TEST_BASE_DIR = os.path.dirname(os.path.dirname(__file__))
    TEST_DATA_DIR = os.path.join(TEST_BASE_DIR, 'tests', 'data')

    def setUp(self):
        self.user = User.objects.create_user('testuser', email='test@example.com', password='test')
        self.user.user_permissions.add(Permission.objects.get(name='Can risk'))
        self.user.save()
        self.token = Token.objects.create(user=self.user)
        self.token.save()
        self.client = APIClient(enforce_csrf_checks=True)
        self.client.credentials(HTTP_AUTHORIZATION='Token ' + self.token.key)

    @pytest.mark.req_WS_CORE_112
    def test_post_to_model_validation_error(self):
        """Test post_to_model method handles ValidationError properly."""
        # Create a mock view that inherits from ModelWebServiceMixin
        class MockView(ModelWebServiceMixin):
            serializer_class = MagicMock()
            any_perms = ['boadicea_auth.can_risk']

        view = MockView()
        request = APIRequestFactory().post('/', {'user_id': 'test', 'pedigree_data': 'test_data'}, format='json')
        request.data = {'user_id': 'test', 'pedigree_data': 'test_data'}
        request.user = self.user

        # Mock serializer to raise ValidationError
        mock_serializer = MagicMock()
        mock_serializer.is_valid.return_value = False
        mock_serializer.errors = {"field": ["error message"]}
        view.serializer_class.return_value = mock_serializer

        response = view.post_to_model(request, settings.BC_MODEL)

        self.assertEqual(response.status_code, status.HTTP_400_BAD_REQUEST)
        self.assertEqual(response.data, {"field": ["error message"]})

    @pytest.mark.req_WS_CORE_112
    @patch('bws.rest_api.Predictions')
    @patch('bws.rest_api.PedigreeFile')
    @patch('bws.rest_api.ModelParams')
    @patch('tempfile.mkdtemp')
    @patch('shutil.rmtree')
    def test_post_to_model_pedigree_error(self, mock_rmtree, mock_mkdtemp, mock_model_params, mock_pedigree_file, mock_predictions):
        """Test post_to_model method handles PedigreeError properly."""
        # Create a mock view
        class MockView(ModelWebServiceMixin):
            serializer_class = MagicMock()
            any_perms = ['boadicea_auth.can_risk']

        view = MockView()
        request = APIRequestFactory().post('/', {'user_id': 'test', 'pedigree_data': 'test_data'}, format='json')
        request.user = self.user

        # Mock valid serializer
        mock_serializer = MagicMock()
        mock_serializer.is_valid.return_value = True
        mock_serializer.validated_data = {'pedigree_data': 'test_data', 'user_id': 'test'}
        view.serializer_class.return_value = mock_serializer

        # Mock PedigreeFile to raise PedigreeError for the first pedigree only
        mock_pf = MagicMock()
        first_pedigree = MagicMock()
        second_pedigree = MagicMock()
        first_pedigree.validateAll.side_effect = PedigreeError("Test pedigree error")
        first_pedigree.get_target.return_value = SimpleNamespace(age='50', pid='1')
        first_pedigree.famid = 'FAM1'
        first_pedigree.hgt = -1
        first_pedigree.mdensity = None
        first_pedigree.ethnicity = None
        first_pedigree.is_ashkn.return_value = False
        first_pedigree.is_carrier_probs_viable.return_value = True

        second_pedigree.validateAll.return_value = []
        second_pedigree.get_target.return_value = SimpleNamespace(age='50', pid='2')
        second_pedigree.famid = 'FAM2'
        second_pedigree.hgt = -1
        second_pedigree.mdensity = None
        second_pedigree.ethnicity = None
        second_pedigree.is_ashkn.return_value = False
        second_pedigree.is_carrier_probs_viable.return_value = True

        mock_pf.pedigrees = [first_pedigree, second_pedigree]
        mock_pedigree_file.return_value = mock_pf

        # Mock Predictions return value for the successful second pedigree
        mock_prediction_result = MagicMock(
            version='1',
            mutation_probabilties=[],
            cancer_risks=[],
            baseline_cancer_risks=[],
            lifetime_cancer_risk=[],
            baseline_lifetime_cancer_risk=[],
            ten_yr_cancer_risk=[],
            baseline_ten_yr_cancer_risk=[],
            ten_yr_nhs_protocol=[]
        )
        mock_predictions.return_value = mock_prediction_result

        # Mock model params
        mock_mp = MagicMock()
        mock_mp.population = 'UK'
        mock_mp.mutation_frequency = {}
        mock_mp.mutation_sensitivity = {}
        mock_mp.cancer_rates = 'UK'
        mock_mp.isashk = False
        mock_model_params.factory.return_value = mock_mp

        # Mock temp directory
        mock_mkdtemp.return_value = '/tmp/test'
        request.data = {'user_id': 'test', 'pedigree_data': mock_pedigree_file}
        response = view.post_to_model(request, settings.BC_MODEL)

        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertIn('errors', response.data)
        self.assertTrue(any('Test pedigree error' in str(err) for err in response.data['errors']))

    @pytest.mark.req_WS_CORE_112
    @patch('bws.rest_api.Predictions')
    @patch('bws.rest_api.PedigreeFile')
    @patch('bws.rest_api.ModelParams')
    @patch('tempfile.mkdtemp')
    @patch('shutil.rmtree')
    def test_post_to_model_model_error(self, mock_rmtree, mock_mkdtemp, mock_model_params, mock_pedigree_file, mock_predictions):
        """Test post_to_model method handles ModelError properly."""
        # Create a mock view
        class MockView(ModelWebServiceMixin):
            serializer_class = MagicMock()
            any_perms = ['boadicea_auth.can_risk']

        view = MockView()
        request = APIRequestFactory().post('/', {'user_id': 'test', 'pedigree_data': 'test_data'}, format='json')
        request.user = self.user

        # Mock valid serializer
        mock_serializer = MagicMock()
        mock_serializer.is_valid.return_value = True
        mock_serializer.validated_data = {'pedigree_data': 'test_data', 'user_id': 'test'}
        view.serializer_class.return_value = mock_serializer

        # Mock PedigreeFile
        mock_pedi = MagicMock()
        mock_pedi.get_target.return_value = SimpleNamespace(age='50', pid='1')
        mock_pedi.famid = 'FAM1'
        mock_pedi.hgt = -1
        mock_pedi.mdensity = None
        mock_pedi.ethnicity = None
        mock_pedi.is_ashkn.return_value = False
        mock_pedi.is_carrier_probs_viable.return_value = True
        mock_pf = MagicMock()
        mock_pf.pedigrees = [mock_pedi]
        mock_pedigree_file.return_value = mock_pf

        # Mock model params
        mock_mp = MagicMock()
        mock_mp.population = 'UK'
        mock_mp.mutation_frequency = {}
        mock_mp.mutation_sensitivity = {}
        mock_mp.cancer_rates = 'UK'
        mock_mp.isashk = False
        mock_model_params.factory.return_value = mock_mp

        # Mock temp directory
        mock_mkdtemp.return_value = '/tmp/test'

        # Mock Predictions to raise ModelError
        mock_predictions.side_effect = ModelError("Test model error")
        request.data = {'user_id': 'test', 'pedigree_data': mock_pedigree_file}
        response = view.post_to_model(request, settings.BC_MODEL)

        self.assertEqual(response.status_code, status.HTTP_400_BAD_REQUEST)
        # Should contain error details
        self.assertIn('Test model error', str(json.loads(response.content)))
