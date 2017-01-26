from django.test import TestCase
from bws.risk_factors import RiskFactors
from boadicea.exceptions import RiskFactorError


class RiskFactorTests(TestCase):

    def test_round_trip(self):
        ''' Test encoding of risk categories and decoding returns same risk categories. '''
        category1 = list(RiskFactors.categories.values())
        category2 = RiskFactors.decode(RiskFactors.encode(category1))
        self.assertListEqual(category1, category2, "round trip encode/decode")

    def test_wrong_no_risks(self):
        ''' Test that an error is raised when the wrong number of risks is submitted '''
        self.assertRaises(RiskFactorError, RiskFactors.encode, [7, 4, 4, 3, 4])
