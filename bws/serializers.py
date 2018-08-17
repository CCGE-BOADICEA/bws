"""
I/O serializers for the web-services.
"""

from rest_framework import serializers
from django.conf import settings
from django.core.files.base import File
from bws.risk_factors import RiskFactors


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


class BaseInputSerializer(serializers.Serializer):
    """ Base serializer for cancer risk calculation input. """
    user_id = serializers.CharField(min_length=4, max_length=40, required=True)
    pedigree_data = FileField()

    @classmethod
    def get_mutation_frequency_field(cls, model):
        """ Get the mutation frequency choice field. """
        return serializers.ChoiceField(choices=list(model['MUTATION_FREQUENCIES'].keys()),
                                       default='UK', help_text="Mutation frequency")

    @classmethod
    def get_gene_mutation_frequency_fields(cls, model):
        """ Get a list of the gene mutation frequency fields strings. """
        MIN_MUT_FREQ = str(settings.MIN_MUTATION_FREQ)
        MAX_MUT_FREQ = str(settings.MAX_MUTATION_FREQ)
        fields = []
        for gene in model['GENES']:
            fields.append(gene.lower() + "_mut_frequency = serializers.FloatField(required=False, "
                          "max_value="+MAX_MUT_FREQ+", min_value="+MIN_MUT_FREQ+")")
        return fields

    @classmethod
    def get_gene_mutation_sensitivity_fields(cls, model):
        """ Get a list of the gene mutation sensitivity fields strings. """
        gts = model['GENETIC_TEST_SENSITIVITY']
        fields = []
        for gene in model['GENES']:
            fields.append(gene.lower() + "_mut_sensitivity = serializers.FloatField(required=False, default=" +
                          str(gts[gene]) + ", max_value=1, min_value=0)")
        return fields

    @classmethod
    def get_cancer_rates_field(cls, model):
        """ Get the cancer rates choice field. """
        return serializers.ChoiceField(choices=list(model['CANCER_RATES'].keys()))


class BwsInputSerializer(BaseInputSerializer):
    """ Boadicea breast cancer input fields. """
    bc_model = settings.BC_MODEL
    mut_freq = BaseInputSerializer.get_mutation_frequency_field(bc_model)

    for f in BaseInputSerializer.get_gene_mutation_frequency_fields(bc_model):
        exec(f)

    for f in BaseInputSerializer.get_gene_mutation_sensitivity_fields(bc_model):
        exec(f)
    cancer_rates = BaseInputSerializer.get_cancer_rates_field(bc_model)


class OwsInputSerializer(BaseInputSerializer):
    """ Boadicea breast cancer input fields. """
    oc_model = settings.OC_MODEL
    mut_freq = BaseInputSerializer.get_mutation_frequency_field(oc_model)

    for f in BaseInputSerializer.get_gene_mutation_frequency_fields(oc_model):
        exec(f)

    for f in BaseInputSerializer.get_gene_mutation_sensitivity_fields(oc_model):
        exec(f)
    cancer_rates = BaseInputSerializer.get_cancer_rates_field(oc_model)


class BwsExtendedInputSerializer(BwsInputSerializer):
    """ Other input parameters. """
    risk_factor_code = serializers.IntegerField(max_value=RiskFactors.get_max_factor(),
                                                min_value=0, default=0)
    prs = serializers.JSONField(required=False)


class OwsExtendedInputSerializer(OwsInputSerializer):
    """ Other input parameters. """
    risk_factor_code = serializers.IntegerField(max_value=RiskFactors.get_max_factor(),
                                                min_value=0, default=0)
    prs = serializers.JSONField(required=False)


class PedigreeResultSerializer(serializers.Serializer):
    family_id = serializers.CharField(read_only=True)
    proband_id = serializers.CharField(read_only=True)
    cancer_risks = serializers.ListField(read_only=True, required=False)
    baseline_cancer_risks = serializers.ListField(read_only=True, required=False)
    lifetime_cancer_risk = serializers.ListField(read_only=True, required=False)
    baseline_lifetime_cancer_risk = serializers.ListField(read_only=True, required=False)
    ten_yr_cancer_risk = serializers.ListField(read_only=True, required=False)
    baseline_ten_yr_cancer_risk = serializers.ListField(read_only=True, required=False)
    mutation_probabilties = serializers.ListField(read_only=True)


class OutputSerializer(serializers.Serializer):
    """ Boadicea result. """
    version = serializers.CharField(read_only=True)
    timestamp = serializers.DateTimeField(read_only=True)
    mutation_frequency = serializers.DictField(read_only=True)
    mutation_sensitivity = serializers.DictField(read_only=True)
    cancer_incidence_rates = serializers.CharField(read_only=True)
    prs = serializers.DictField(read_only=True, required=False)
    risk_factors = serializers.DictField(read_only=True, required=False)
    warnings = serializers.ListField(read_only=True, required=False)
    pedigree_result = PedigreeResultSerializer(read_only=True, many=True)
