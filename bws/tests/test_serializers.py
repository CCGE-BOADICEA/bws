"""
Unit tests for bws.serializers — field-level validation and edge cases.

© 2026 University of Cambridge
SPDX-FileCopyrightText: 2026 University of Cambridge
SPDX-License-Identifier: GPL-3.0-or-later
"""
import io
import pytest

from django.conf import settings
from django.core.files.base import ContentFile
from django.test import TestCase
from rest_framework import serializers

from bws.serializers import (
    BaseInputSerializer,
    BwsInputSerializer,
    FileField,
    OwsInputSerializer,
    PRSField,
    PwsInputSerializer,
)


class FileFieldTests(TestCase):
    """ Tests for FileField.to_internal_value(). """

    def setUp(self):
        self.field = FileField()

    @pytest.mark.req_WS_CORE_113
    def test_string_input_returned_unchanged(self):
        """ A plain string is returned as-is. """
        data = "##CanRisk 1.0\nsome pedigree data"
        result = self.field.to_internal_value(data)
        self.assertEqual(result, data)

    @pytest.mark.req_WS_CORE_113
    def test_file_object_returns_content_as_string(self):
        """ An uploaded file object is read and returned as a string. """
        content = b"##CanRisk 1.0\nsome pedigree data"
        uploaded = ContentFile(content, name="pedigree.canrisk")
        result = self.field.to_internal_value(uploaded)
        self.assertEqual(result, content.decode())


class BaseInputSerializerExtraFieldTests(TestCase):
    """ Tests for BaseInputSerializer.is_valid() extra-field rejection. """

    @pytest.mark.req_WS_CORE_114
    def test_extra_field_raises_validation_error(self):
        """ Submitting a field not in the serializer raises ValidationError. """
        data = {
            'user_id': 'testuser',
            'pedigree_data': '##CanRisk 1.0',
            'unknown_extra_field': 'bad_value',
        }
        ser = BwsInputSerializer(data=data)
        with self.assertRaises(serializers.ValidationError) as ctx:
            ser.is_valid(raise_exception=True)
        self.assertIn('Input Field Error', str(ctx.exception.detail))

    @pytest.mark.req_WS_CORE_114
    def test_no_extra_fields_does_not_raise(self):
        """ Submitting only recognised fields does not raise a ValidationError for extra fields. """
        data = {
            'user_id': 'testuser',
            'pedigree_data': '##CanRisk 1.0',
            'mut_freq': 'UK',
            'cancer_rates': 'UK',
        }
        ser = BwsInputSerializer(data=data)
        try:
            ser.is_valid(raise_exception=True)
        except serializers.ValidationError as exc:
            if 'Input Field Error' in str(exc.detail):
                self.fail("is_valid raised Input Field Error for a known field")


class BwsInputSerializerFieldTests(TestCase):
    """ Tests that BwsInputSerializer exposes the correct model-specific fields. """

    @pytest.mark.req_WS_CORE_115
    def test_cancer_rates_field_present(self):
        """ cancer_rates field exists in BwsInputSerializer. """
        ser = BwsInputSerializer()
        self.assertIn('cancer_rates', ser.fields)

    @pytest.mark.req_WS_CORE_115
    def test_mut_freq_field_present(self):
        """ mut_freq field exists in BwsInputSerializer. """
        ser = BwsInputSerializer()
        self.assertIn('mut_freq', ser.fields)

    @pytest.mark.req_WS_CORE_115
    def test_gene_sensitivity_fields_present(self):
        """ A sensitivity field exists for every BC model gene. """
        ser = BwsInputSerializer()
        for gene in settings.BC_MODEL['GENES']:
            self.assertIn(gene.lower() + '_mut_sensitivity', ser.fields)

    @pytest.mark.req_WS_CORE_115
    def test_uk_is_valid_cancer_rates_choice(self):
        """ 'UK' is a valid cancer_rates choice. """
        ser = BwsInputSerializer()
        choices = list(ser.fields['cancer_rates'].choices.keys())
        self.assertIn('UK', choices)


class OwsInputSerializerFieldTests(TestCase):
    """ Tests that OwsInputSerializer exposes OC model fields. """

    @pytest.mark.req_WS_CORE_115
    def test_oc_gene_sensitivity_fields_present(self):
        """ A sensitivity field exists for every OC model gene. """
        ser = OwsInputSerializer()
        for gene in settings.OC_MODEL['GENES']:
            self.assertIn(gene.lower() + '_mut_sensitivity', ser.fields)


class PwsInputSerializerFieldTests(TestCase):
    """ Tests that PwsInputSerializer exposes PC model fields. """

    @pytest.mark.req_WS_CORE_115
    def test_pc_gene_sensitivity_fields_present(self):
        """ A sensitivity field exists for every PC model gene. """
        ser = PwsInputSerializer()
        for gene in settings.PC_MODEL['GENES']:
            self.assertIn(gene.lower() + '_mut_sensitivity', ser.fields)


class GetMutationFrequencyFieldTests(TestCase):
    """ Tests for BaseInputSerializer.get_mutation_frequency_field(). """

    @pytest.mark.req_WS_CORE_113
    def test_returns_choice_field_with_model_keys(self):
        """ Returns a ChoiceField whose choices match the model's MUTATION_FREQUENCIES keys. """
        model = settings.BC_MODEL
        field = BaseInputSerializer.get_mutation_frequency_field(model)
        self.assertIsInstance(field, serializers.ChoiceField)
        expected_choices = list(model['MUTATION_FREQUENCIES'].keys())
        actual_choices = list(field.choices.keys())
        for key in expected_choices:
            self.assertIn(key, actual_choices)
