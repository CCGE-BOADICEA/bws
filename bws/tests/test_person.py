"""
Unit tests for bws.person — Person, Male, Female classes.

Focuses on logic not already covered by tests_pedigree_validation.py:
is_complete(), is_target(), sex(), constructor behaviour, and validation
of the dead and ashkn flags.

© 2026 University of Cambridge
SPDX-FileCopyrightText: 2026 University of Cambridge
SPDX-License-Identifier: GPL-3.0-or-later
"""
from copy import deepcopy
import os
import pytest

from bws.cancer import Cancers, Cancer
from bws.exceptions import PersonError
from bws.pedigree_file import PedigreeFile
from bws.person import Male, Female, Person
from django.conf import settings
from django.test import TestCase


TEST_BASE_DIR = os.path.dirname(os.path.dirname(__file__))
TEST_DATA_DIR = os.path.join(TEST_BASE_DIR, 'tests', 'data')


def _load_pedigree():
    """ Return a fresh PedigreeFile and its first pedigree from d3.bwa. """
    with open(os.path.join(TEST_DATA_DIR, 'd3.bwa'), 'r') as f:
        pedigree_file = PedigreeFile(f.read())
    return pedigree_file


class PersonSexTests(TestCase):
    """ Tests for Male.sex() and Female.sex(). """

    @pytest.mark.req_WS_VALIDATION_255
    def test_male_sex(self):
        """ Male.sex() returns 'M'. """
        m = Male('FAM1', 'M1', 'P01', '0', '0')
        self.assertEqual(m.sex(), 'M')

    @pytest.mark.req_WS_VALIDATION_255
    def test_female_sex(self):
        """ Female.sex() returns 'F'. """
        f = Female('FAM1', 'F1', 'P02', '0', '0')
        self.assertEqual(f.sex(), 'F')


class PersonConstructorTests(TestCase):
    """ Tests for Person.__init__ constructor edge cases. """

    @pytest.mark.req_WS_VALIDATION_256
    def test_famid_hyphen_removed(self):
        """ Hyphens in famid are stripped on construction. """
        p = Male('FAM-1', 'M1', 'P01', '0', '0')
        self.assertEqual(p.famid, 'FAM1')

    @pytest.mark.req_WS_VALIDATION_256
    def test_famid_truncated_to_8_chars(self):
        """ famid is truncated to 8 characters. """
        p = Male('ABCDEFGHIJ', 'M1', 'P01', '0', '0')
        self.assertEqual(len(p.famid), 8)

    @pytest.mark.req_WS_VALIDATION_256
    def test_name_truncated_to_8_chars(self):
        """ name is truncated to 8 characters. """
        p = Male('FAM1', 'VERYLONGNAME', 'P01', '0', '0')
        self.assertEqual(len(p.name), 8)


class PersonIsCompleteTests(TestCase):
    """ Tests for Person.is_complete(). """

    @pytest.mark.req_WS_VALIDATION_253
    def test_missing_age_and_yob_is_incomplete(self):
        """ Returns False when both age and yob are '0'. """
        p = Female('FAM1', 'F1', 'P01', '0', '0', age='0', yob='0')
        self.assertFalse(p.is_complete())

    @pytest.mark.req_WS_VALIDATION_253
    def test_missing_yob_only_is_incomplete(self):
        """ Returns False when yob is '0' regardless of age having a value. """
        p = Female('FAM1', 'F1', 'P01', '0', '0', age='50', yob='0')
        self.assertFalse(p.is_complete())

    @pytest.mark.req_WS_VALIDATION_253
    def test_missing_age_only_without_cancer_is_incomplete(self):
        """ Returns False when age is '0' regardless of yob having a value. """
        p = Female('FAM1', 'F1', 'P01', '0', '0', age='0', yob='1970')
        self.assertFalse(p.is_complete())

    @pytest.mark.req_WS_VALIDATION_253
    def test_complete_person_with_age_and_yob(self):
        """ Returns True when both age and yob are provided. """
        p = Female('FAM1', 'F1', 'P01', '0', '0', age='50', yob='1970')
        self.assertTrue(p.is_complete())

    @pytest.mark.req_WS_VALIDATION_253
    def test_incomplete_but_has_cancer_is_complete(self):
        """ Returns True when age/yob unknown but the person has a cancer diagnosis. """
        cancers = Cancers(bc1=Cancer('40'))
        p = Female('FAM1', 'F1', 'P01', '0', '0', age='0', yob='0', cancers=cancers)
        self.assertTrue(p.is_complete())


class PersonIsTargetTests(TestCase):
    """ Tests for Person.is_target(). """

    @pytest.mark.req_WS_VALIDATION_254
    def test_target_flag_zero_is_not_target(self):
        """ Returns False when target is '0'. """
        p = Female('FAM1', 'F1', 'P01', '0', '0', target='0')
        self.assertFalse(p.is_target())

    @pytest.mark.req_WS_VALIDATION_254
    def test_target_flag_one_is_target(self):
        """ Returns True when target is '1'. """
        p = Female('FAM1', 'F1', 'P01', '0', '0', target='1')
        self.assertTrue(p.is_target())


class PersonDeadFlagValidationTests(TestCase):
    """ Tests that the dead flag is validated by Person.validate(). """

    def setUp(self):
        self.pedigree_file = _load_pedigree()

    @pytest.mark.req_WS_VALIDATION_257
    def test_invalid_dead_flag_raises_person_error(self):
        """ A dead value other than '0' or '1' raises PersonError. """
        pf = deepcopy(self.pedigree_file)
        apedigree = pf.pedigrees[0]
        person = apedigree.get_person_by_name('F1')
        person.dead = '2'
        with self.assertRaisesRegex(PersonError, r"invalid vital status"):
            person.validate(apedigree)

    @pytest.mark.req_WS_VALIDATION_257
    def test_dead_zero_is_valid(self):
        """ dead='0' (alive) passes validation. """
        pf = deepcopy(self.pedigree_file)
        apedigree = pf.pedigrees[0]
        person = apedigree.get_person_by_name('F1')
        person.dead = '0'
        person.validate(apedigree)

    @pytest.mark.req_WS_VALIDATION_257
    def test_dead_one_is_valid(self):
        """ dead='1' (deceased) passes validation. """
        pf = deepcopy(self.pedigree_file)
        apedigree = pf.pedigrees[0]
        person = apedigree.get_person_by_name('F1')
        person.dead = '1'
        person.validate(apedigree)


class PersonAshknFlagValidationTests(TestCase):
    """ Tests that the Ashkenazi flag is validated by Person.validate(). """

    def setUp(self):
        self.pedigree_file = _load_pedigree()

    @pytest.mark.req_WS_VALIDATION_258
    def test_invalid_ashkn_flag_raises_person_error(self):
        """ An ashkn value other than '0' or '1' raises PersonError. """
        pf = deepcopy(self.pedigree_file)
        apedigree = pf.pedigrees[0]
        person = apedigree.get_person_by_name('F1')
        person.ashkn = '2'
        with self.assertRaisesRegex(PersonError, r"invalid Ashkenazi origin parameter"):
            person.validate(apedigree)

    @pytest.mark.req_WS_VALIDATION_258
    def test_ashkn_zero_is_valid(self):
        """ ashkn='0' passes validation. """
        pf = deepcopy(self.pedigree_file)
        apedigree = pf.pedigrees[0]
        person = apedigree.get_person_by_name('F1')
        person.ashkn = '0'
        person.validate(apedigree)

    @pytest.mark.req_WS_VALIDATION_258
    def test_ashkn_one_is_valid(self):
        """ ashkn='1' passes validation. """
        pf = deepcopy(self.pedigree_file)
        apedigree = pf.pedigrees[0]
        person = apedigree.get_person_by_name('F1')
        person.ashkn = '1'
        person.validate(apedigree)