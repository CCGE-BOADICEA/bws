"""
BOADICEA web-service testing.

© 2026 University of Cambridge
SPDX-FileCopyrightText: 2026 University of Cambridge
SPDX-License-Identifier: GPL-3.0-or-later

Unit tests for bws.cancer module.
"""
from unittest.mock import MagicMock, patch

from django.test import TestCase

from bws.cancer import (
    Cancer,
    Genes,
    PathologyTest,
    PathologyError,
    GeneticTest,
    Cancers,
    CancerError,
)


class DummyPerson:
    def __init__(self, pid='id1', famid='fam', age='50', yob='1970', sex='F'):
        self.pid = pid
        self.famid = famid
        self.age = age
        self.yob = yob
        self._sex = sex
        self.pathology = PathologyTest.factory_default()
        self.gtests = []
        self.cancers = Cancers()
    def sex(self):
        return self._sex


class TestCancerUnit(TestCase):

    @patch('bws.cancer.settings')
    def test_genes_unique_and_all(self, mock_settings):
        """Test that unique and all genes are correctly identified from settings."""
        mock_settings.BC_MODEL = {'GENES': ['A', 'B', 'C']}
        mock_settings.OC_MODEL = {'GENES': ['C', 'D', 'E']}
        mock_settings.PC_MODEL = {'GENES': ['F', 'B', 'G']}

        oc_unique = Genes.get_unique_oc_genes()
        pc_unique = Genes.get_unique_pc_genes()
        all_genes = Genes.get_all_model_genes()

        self.assertEqual(set(oc_unique), {'D', 'E'})
        self.assertEqual(set(pc_unique), {'F', 'G'})
        self.assertIsInstance(all_genes, list)

    def test_pathology_write_default_and_er_positive(self):
        """Test that PathologyTest.write correctly formats default and ER positive cases."""
        # default all unknown -> -1
        tests = PathologyTest.factory_default()
        self.assertEqual(PathologyTest.write(tests).strip(), '-1')

        # ER positive -> code 0
        tests = tests._replace(er=PathologyTest(PathologyTest.ESTROGEN_RECEPTOR_TEST, 'P'))
        self.assertEqual(PathologyTest.write(tests).strip(), '0')

    def test_pathology_validate_invalid_status_raises(self):
        """Test that invalid pathology status raises PathologyError."""
        person = DummyPerson()
        person.pathology = person.pathology._replace(er=PathologyTest('1', 'X'))
        # actual message: "has been assigned an invalid ... status"
        with self.assertRaisesRegex(PathologyError, "invalid"):
            PathologyTest.validate(person)

    def test_pathology_validate_requires_bc_for_result(self):
        """Test that a pathology result without a BC diagnosis raises PathologyError."""
        person = DummyPerson()
        person.pathology = person.pathology._replace(er=PathologyTest('1', 'P'))
        # person.cancers defaults to Cancers() i.e. all unaffected — no bc1 needed
        # actual message: "has not developed breast cancer but has been assigned a breast cancer pathology test result"
        with self.assertRaisesRegex(PathologyError, "not developed breast cancer"):
            PathologyTest.validate(person)

    def test_genetictest_get_data_and_compare(self):
        """Test that GeneticTest.get_genetic_test_data returns correct codes and compareTestResults works."""
        self.assertEqual(GeneticTest('S', 'N').get_genetic_test_data(), '0')
        self.assertEqual(GeneticTest('S', 'P').get_genetic_test_data(), '1')
        self.assertEqual(GeneticTest('T', 'N').get_genetic_test_data(), '2')
        self.assertEqual(GeneticTest('T', 'P').get_genetic_test_data(), '3')

        # HOXB13 behaviour
        self.assertEqual(GeneticTest('0', 'HET', isHOXB13=True).get_genetic_test_data(), '6')
        self.assertEqual(GeneticTest('0', 'HOM', isHOXB13=True).get_genetic_test_data(), '7')

        # untested
        self.assertEqual(GeneticTest('0', '0',   isHOXB13=True).get_genetic_test_data(), '-1')

        # compareTestResults
        p1, p2 = MagicMock(), MagicMock()
        p1.gtests = [GeneticTest('S', 'N')]
        p2.gtests = [GeneticTest('S', 'N')]
        self.assertTrue(GeneticTest.compareTestResults(p1, p2))
        p2.gtests = [GeneticTest('S', 'P')]
        self.assertFalse(GeneticTest.compareTestResults(p1, p2))

    @patch('bws.cancer.settings')
    def test_cancers_validate_age_and_sex_rules(self, mock_settings):
        """Test that Cancers.validate raises CancerError for invalid age formats,
        ages above MAX_AGE, and sex-specific cancer rules."""
        mock_settings.MAX_AGE = 125

        # invalid age format
        person = DummyPerson()
        person.cancers = Cancers(bc1=Cancer('abc'))
        with self.assertRaisesRegex(CancerError, "age at cancer diagnosis"):
            Cancers.validate(person)

        # age at diagnosis > MAX_AGE
        person = DummyPerson(age='50')
        person.cancers = Cancers(bc1=Cancer('200'))
        with self.assertRaisesRegex(CancerError, "specified with an integer in the range 1"):
            Cancers.validate(person)

        # male with ovarian cancer
        person = DummyPerson(sex='M')
        person.cancers = Cancers(oc=Cancer('45'))
        with self.assertRaisesRegex(CancerError, "is male but has been assigned an ovarian cancer diagnosis"):
            Cancers.validate(person)

        # female with prostate cancer
        person = DummyPerson(sex='F')
        person.cancers = Cancers(prc=Cancer('30'))
        with self.assertRaisesRegex(CancerError, "female but has been assigned an prostate cancer diagnosis"):
            Cancers.validate(person)
