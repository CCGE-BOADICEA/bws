""" Command line utility. """
from django.core.management.base import BaseCommand
from bws.risk_factors import RiskFactors


class Command(BaseCommand):
    help = 'Decode risk factors into categories, e.g ./manage.py decode_risk_factors 112005'

    def add_arguments(self, parser):
        parser.add_argument('factor', type=int)

    def handle(self, *args, **options):
        risk_factor_code = options['factor']
        categories = RiskFactors.decode(risk_factor_code)
        for idx, cat in enumerate(categories):
            name = RiskFactors.risk_factors[idx].__name__
            print(name + " idx: " + str(cat) + " category: " + RiskFactors.risk_factors[idx].cats[cat])
