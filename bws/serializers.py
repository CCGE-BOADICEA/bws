"""
I/O serializers for the web-services.

© 2023 University of Cambridge
SPDX-FileCopyrightText: 2023 University of Cambridge
SPDX-License-Identifier: GPL-3.0-or-later
"""

from rest_framework import serializers
from django.conf import settings
from django.core.files.base import File
from rest_framework.fields import JSONField


class FileField(serializers.Field):
    """
    Pedigree field object serialized into a string representation. The field can
    be a str or and uploaded file type.
    """
    def to_representation(self, obj):
        return obj

    def to_internal_value(self, obj):
        assert(isinstance(obj, str) or isinstance(obj, File))

        if isinstance(obj, File):
            pedigree_data = ''
            for chunk in obj.chunks():
                pedigree_data += chunk.decode("utf-8")
            return pedigree_data
        else:
            return obj


class PRSField(JSONField):
    alpha = serializers.FloatField(required=True, min_value=0, max_value=1)
    zscore = serializers.FloatField(required=True)
    percent = serializers.FloatField(required=False)


class BaseInputSerializer(serializers.Serializer):
    """ Base serializer for cancer risk calculation input. """
    user_id = serializers.CharField(min_length=4, max_length=40, required=True, help_text="Unique end user ID")
    pedigree_data = FileField(help_text="CanRisk pedigree data file")
    prs = PRSField(label="Polygenic Risk Score",
                   help_text='JSON of alpha and zscore values, e.g. {"alpha":0.501,"zscore":1.65}', required=False)

    @classmethod
    def get_mutation_frequency_field(cls, model):
        """ Get the mutation frequency choice field. """
        return serializers.ChoiceField(choices=list(model['MUTATION_FREQUENCIES'].keys()),
                                       default='UK', help_text="Mutation frequency")

    # @classmethod
    # def get_gene_mutation_frequency_fields(cls, model):
    #     """ Get a list of the gene mutation frequency fields strings. """
    #     MIN_MUT_FREQ = str(settings.MIN_MUTATION_FREQ)
    #     MAX_MUT_FREQ = str(settings.MAX_MUTATION_FREQ)
    #     fields = []
    #     for gene in model['GENES']:
    #         fields.append(gene.lower() + "_mut_frequency = serializers.FloatField(required=False, "
    #                       "max_value="+MAX_MUT_FREQ+", min_value="+MIN_MUT_FREQ+")")
    #     return fields

    @classmethod
    def get_gene_mutation_sensitivity_fields(cls, model):
        """ Get a list of the gene mutation sensitivity fields strings. """
        gts = model['GENETIC_TEST_SENSITIVITY']
        fields = []
        for gene in model['GENES']:
            fields.append(gene.lower() + "_mut_sensitivity = serializers.FloatField(required=False, default=" +
                          str(gts[gene]) + ", max_value=1, min_value=0, help_text='"+gene+" test sensitivity')")
        return fields

    @classmethod
    def get_cancer_rates_field(cls, model):
        """ Get the cancer rates choice field. """
        return serializers.ChoiceField(choices=list(model['CANCER_RATES'].keys()), help_text="Cancer incidence rates")

    def is_valid(self, raise_exception=True):
        """ Check that superfluous input flags have not been included. """
        if hasattr(self, 'initial_data'):
            inp_keys = self.initial_data.keys()     # all the user input keys
            ser_keys = self.fields.keys()           # all the serializer fields
            extra_fields = list(filter(lambda key: key not in ser_keys, inp_keys))
            if len(extra_fields) > 0:
                msg = ", ".join(extra_fields)
                raise serializers.ValidationError({'Input Field Error': f'Extra input field(s) found: {msg}'})
        return super(BaseInputSerializer, self).is_valid(raise_exception=raise_exception)


class BwsInputSerializer(BaseInputSerializer):
    """ Boadicea breast cancer input fields. """
    bc_model = settings.BC_MODEL
    mut_freq = BaseInputSerializer.get_mutation_frequency_field(bc_model)

    # for f in BaseInputSerializer.get_gene_mutation_frequency_fields(bc_model):
    #     exec(f)

    for f in BaseInputSerializer.get_gene_mutation_sensitivity_fields(bc_model):
        exec(f)
    cancer_rates = BaseInputSerializer.get_cancer_rates_field(bc_model)


class OwsInputSerializer(BaseInputSerializer):
    """ Ovarian cancer input fields. """
    oc_model = settings.OC_MODEL
    mut_freq = BaseInputSerializer.get_mutation_frequency_field(oc_model)

    # for f in BaseInputSerializer.get_gene_mutation_frequency_fields(oc_model):
    #     exec(f)

    for f in BaseInputSerializer.get_gene_mutation_sensitivity_fields(oc_model):
        exec(f)
    cancer_rates = BaseInputSerializer.get_cancer_rates_field(oc_model)


class PwsInputSerializer(BaseInputSerializer):
    """ Prostate cancer input fields. """
    pc_model = settings.PC_MODEL
    mut_freq = BaseInputSerializer.get_mutation_frequency_field(pc_model)

    for f in BaseInputSerializer.get_gene_mutation_sensitivity_fields(pc_model):
        exec(f)
    cancer_rates = BaseInputSerializer.get_cancer_rates_field(pc_model)


class PedigreeResultSerializer(serializers.Serializer):
    """ Cancer model risks and mutation probabilities for a pedigree. """
    family_id = serializers.CharField(read_only=True, help_text="Family ID")
    proband_id = serializers.CharField(read_only=True, help_text="Proband ID")
    ethnicity = serializers.CharField(read_only=True, help_text="BioBank ethnicity (UK only)")
    ons_ethnicity = serializers.CharField(read_only=True, help_text="ONS ethnicity (UK only)")
    mutation_frequency = serializers.DictField(read_only=True)
    risk_factors = serializers.DictField(read_only=True, required=False)
    prs = serializers.DictField(read_only=True, required=False)
    cancer_risks = serializers.ListField(read_only=True, required=False)
    baseline_cancer_risks = serializers.ListField(read_only=True, required=False)
    lifetime_cancer_risk = serializers.ListField(read_only=True, required=False)
    baseline_lifetime_cancer_risk = serializers.ListField(read_only=True, required=False)
    ten_yr_cancer_risk = serializers.ListField(read_only=True, required=False)
    ten_yr_nhs_protocol = serializers.ListField(read_only=True, required=False)
    baseline_ten_yr_cancer_risk = serializers.ListField(read_only=True, required=False)
    mutation_probabilties = serializers.ListField(read_only=True)


class OutputSerializer(serializers.Serializer):
    """ Cancer model results. """
    version = serializers.CharField(read_only=True)
    timestamp = serializers.DateTimeField(read_only=True)
    mutation_frequency = serializers.DictField(read_only=True)
    mutation_sensitivity = serializers.DictField(read_only=True)
    cancer_incidence_rates = serializers.CharField(read_only=True)
    warnings = serializers.ListField(read_only=True, required=False)
    pedigree_result = PedigreeResultSerializer(read_only=True, many=True)
    nocalc = serializers.CharField(read_only=True)


class CombinedInputSerializer(serializers.Serializer):
    ''' Results from ovarian, breast and prostate cancer models. '''
    ows_result = serializers.JSONField(required=True)
    bws_result = serializers.JSONField(required=True)
    pws_result = serializers.JSONField(required=False)


class CombinedOutputSerializer(serializers.Serializer):
    """ Results from ovarian, breast and prostate cancer models. """
    ows_result = OutputSerializer(read_only=True)
    bws_result = OutputSerializer(read_only=True)
    pws_result = OutputSerializer(read_only=True, required=False)
