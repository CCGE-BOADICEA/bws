"""
Command line utility.

© 2023 University of Cambridge
SPDX-FileCopyrightText: 2023 University of Cambridge
SPDX-License-Identifier: GPL-3.0-or-later
"""
from django.core.management.base import BaseCommand
from bws.risk_factors.bc import BCRiskFactors


class Command(BaseCommand):
    help = 'Decode risk factors into categories, e.g ./manage.py decode_risk_factors 112005'

    def add_arguments(self, parser):
        parser.add_argument('factor', type=int)

    def handle(self, *args, **options):
        risk_factor_code = options['factor']
        categories = BCRiskFactors.decode(risk_factor_code)
        for idx, cat in enumerate(categories):
            name = BCRiskFactors.risk_factors[idx].__name__
            print(name + " idx: " + str(cat) + " category: " + BCRiskFactors.risk_factors[idx].cats[cat])
