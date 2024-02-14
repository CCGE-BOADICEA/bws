'''
Prostate cancer risk factors.

© 2024 University of Cambridge
SPDX-FileCopyrightText: 2023 University of Cambridge
SPDX-License-Identifier: GPL-3.0-or-later
'''

from bws.risk_factors.rfs import RiskFactors
from collections import OrderedDict


class PCRiskFactors(RiskFactors):
    ''' Each risk factor for an individual is defined in terms of a category they are in.
        If a factor is unobserved, missing or not applicable, it is assigned category 0,
        and is not taken into account in the calculation. Otherwise a non-zero number is given
        depending on which group they belong to. These are then combined into a single
        risk factor code (see encode() function) that is used by the BOADICEA risk calculation. '''

    # list of risk factor classes
    risk_factors = []

    # dictionary of risk factor name and number of categories
    categories = OrderedDict(
        (rf.snake_name(), len(rf.cats)-1) for rf in risk_factors
    )
