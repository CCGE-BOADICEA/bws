"""
Additional edge case tests for pedigree.py to improve coverage.

© 2026 University of Cambridge
SPDX-FileCopyrightText: 2026 University of Cambridge
SPDX-License-Identifier: GPL-3.0-or-later
"""
import os
import tempfile
import shutil
import pytest
from copy import deepcopy
from datetime import date

from django.test import TestCase, override_settings
from django.conf import settings

from bws.cancer import Cancers, Cancer, BWSGeneticTests, CanRiskGeneticTests, GeneticTest
from bws.pedigree import BwaPedigree, CanRiskPedigree
from bws.pedigree_file import PedigreeFile
from bws.person import Male, Female
from bws.exceptions import PedigreeError
from bws.risk_factors.mdensity import Birads


class PedigreeEdgeCasesTests(TestCase):
    """Test edge cases in pedigree methods."""

    def setUp(self):
        """Set up test fixtures."""
        self.year = date.today().year
        target = Female("FAM1", "F0", "001", "002", "003", target="1", age="30",
                       yob=str(self.year-30), cancers=Cancers())
        self.pedigree = BwaPedigree(people=[target])
        self.pedigree.add_parents(target)

    def test_get_person_by_name_not_found(self):
        """Test getting person by name when not found."""
        result = self.pedigree.get_person_by_name("NonExistent")
        self.assertIsNone(result)

    def test_get_person_not_found(self):
        """Test getting person by ID when not found."""
        result = self.pedigree.get_person("999")
        self.assertIsNone(result)

    def test_get_siblings_no_parents(self):
        """Test getting siblings when person has no parents specified."""
        orphan = Female("FAM1", "F2", "004", "0", "0", target="0", age="25",
                       yob=str(self.year-25), cancers=Cancers())
        self.pedigree.people.append(orphan)
        
        siblings, siblings_yob = self.pedigree.get_siblings(orphan)
        self.assertEqual(len(siblings), 0)
        self.assertEqual(len(siblings_yob), 0)

    def test_write_pedigree_file_birads(self):
        """Test writing input pedigree file for fortran when mdensity has '4' value."""
        cwd = tempfile.mkdtemp(prefix="TEST_", dir="/tmp")
        try:
            filepath = os.path.join(cwd, "test_birads4.ped")
            
            # Create mdensity with NA value
            birads_na = Birads("4")
            result = self.pedigree.write_pedigree_file(
                filepath=filepath,
                mdensity=birads_na,
                model_settings=settings.BC_MODEL
            )
            
            self.assertTrue(os.path.exists(result))
            # Check file contents
            with open(result, 'r') as f:
                content = f.read()
                self.assertIn("00000004", content)
        finally:
            shutil.rmtree(cwd)


class PedigreeInitialization(TestCase):
    """Test edge cases in pedigree initialization."""

    def setUp(self):
        """Set up test fixtures."""
        self.year = date.today().year
        self.cwd = tempfile.mkdtemp(prefix="TEST_", dir="/tmp")
        
    def tearDown(self):
        """Clean up test fixtures."""
        TestCase.tearDown(self)
        shutil.rmtree(self.cwd)

    def test_pedigree_from_records_basic(self):
        """Test creating pedigree from pedigree records."""
        pedigree_data = """FAM001\tpa\t0\tm21\t0\t0\tM\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0:0:0:0\n
FAM001\tma\t0\tf21\t0\t0\tF\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0:0:0:0\n
FAM001\tme\t1\tch1\tm21\tf21\tF\t0\t0\t45\t1980\t0\t0\t0\t0\t0\t0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0:0:0:0"""

        records = [line for line in pedigree_data.strip().split('\n') if line]
        pedigree = CanRiskPedigree(pedigree_records=records, file_type="canrisk4", delim="\t")
        
        self.assertEqual(pedigree.famid, "FAM001")
        self.assertEqual(len(pedigree.people), 3)

    def test_pedigree_invalid_target_value(self):
        """Test that invalid target value raises error."""
        pedigree_data = """FAM001\tpa\t0\tm21\t0\t0\tM\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0:0:0:0\n
FAM001\tma\t0\tf21\t0\t0\tF\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0:0:0:0\n
FAM001\tme\t2\tch1\tm21\tf21\tF\t0\t0\t45\t1980\t0\t0\t0\t0\t0\t0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0:0:0:0"""

        records = [line for line in pedigree_data.strip().split('\n') if line]
        with self.assertRaisesRegex(PedigreeError, r"Target.*'2'.*must be set to"):
            pedigree = CanRiskPedigree(pedigree_records=records, file_type="canrisk4", delim="\t")

    def test_pedigree_duplicate_person_id(self):
        """Test that duplicate person IDs raise error."""
        pedigree_data = """FAM001\tpa\t0\tch1\t0\t0\tM\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0:0:0:0\n
FAM001\tma\t0\tf21\t0\t0\tF\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0:0:0:0\n
FAM001\tme\t1\tch1\tm21\tf21\tF\t0\t0\t45\t1980\t0\t0\t0\t0\t0\t0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0:0:0:0"""

        records = [line for line in pedigree_data.strip().split('\n') if line]
        with self.assertRaisesRegex(PedigreeError, r"Individual ID 'ch1' appears more than once"):
            pedigree = CanRiskPedigree(pedigree_records=records, file_type="canrisk4", delim="\t")

    def test_pedigree_no_target(self):
        """Test that pedigree with no target raises error."""
        pedigree_data = """FAM001\tpa\t0\tm21\t0\t0\tM\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0:0:0:0\n
FAM001\tma\t0\tf21\t0\t0\tF\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0:0:0:0\n
FAM001\tme\t0\tch1\tm21\tf21\tF\t0\t0\t45\t1980\t0\t0\t0\t0\t0\t0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0:0:0:0"""
  
        records = [line for line in pedigree_data.strip().split('\n') if line]
        with self.assertRaisesRegex(PedigreeError, r"has either no index or more than 1 index"):
            pedigree = CanRiskPedigree(pedigree_records=records, file_type="canrisk4", delim="\t")

    def test_pedigree_multiple_targets(self):
        """Test that pedigree with multiple targets raises error."""
        pedigree_data = """FAM001\tpa\t0\tm21\t0\t0\tM\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0:0:0:0\n
FAM001\tma\t1\tf21\t0\t0\tF\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0:0:0:0\n
FAM001\tme\t1\tch1\tm21\tf21\tF\t0\t0\t45\t1980\t0\t0\t0\t0\t0\t0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0:0:0:0"""
  
        records = [line for line in pedigree_data.strip().split('\n') if line]
        with self.assertRaisesRegex(PedigreeError, r"has either no index or more than 1 index"):
            pedigree = CanRiskPedigree(pedigree_records=records, file_type="canrisk4", delim="\t")

    def test_pedigree_with_people_and_records(self):
        """Test creating pedigree with both people and records."""
        # Create a pedigree from records
        pedigree_data = """FAM001\tpa\t0\tm21\t0\t0\tM\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0:0:0:0\n
FAM001\tma\t0\tf21\t0\t0\tF\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0:0:0:0\n
FAM001\tme\t1\tch1\tm21\tf21\tF\t0\t0\t45\t1980\t0\t0\t0\t0\t0\t0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0:0:0:0"""

        records = [line for line in pedigree_data.strip().split('\n') if line]
        
        # Create additional people
        extra_person = Male("FAM001", "M3", "004", "0", "0", age="40", 
                           yob="1960", cancers=Cancers())
        
        # Create pedigree with both
        pedigree = CanRiskPedigree(pedigree_records=records, file_type="canrisk4", delim="\t", people=[extra_person])
        
        self.assertEqual(len(pedigree.people), 4)
        self.assertEqual(pedigree.famid, "FAM001")


class PedigreeValidation(TestCase):
    """Test edge cases in pedigree validation."""

    def setUp(self):
        """Set up test fixtures."""
        self.year = date.today().year

    def test_is_ashkn_with_ashkenazi(self):
        """Test detecting Ashkenazi ancestry."""
        target = Female("FAM1", "F0", "001", "0", "0", target="1", age="30",
                       yob=str(self.year-30), cancers=Cancers(), ashkn="1")
        other = Male("FAM1", "M0", "002", "0", "0", age="60", yob="1950", 
                    cancers=Cancers(), ashkn="0")
        
        pedigree = BwaPedigree(people=[target, other])
        self.assertTrue(pedigree.is_ashkn())

    def test_is_ashkn_no_ashkenazi(self):
        """Test detecting no Ashkenazi ancestry."""
        target = Female("FAM1", "F0", "001", "0", "0", target="1", age="30",
                       yob=str(self.year-30), cancers=Cancers(), ashkn="0")
        other = Male("FAM1", "M0", "002", "0", "0", age="60", yob="1950",
                    cancers=Cancers(), ashkn="0")
        
        pedigree = BwaPedigree(people=[target, other])
        self.assertFalse(pedigree.is_ashkn())


class WritePedigreeFile(TestCase):
    """Additional tests for write methods."""

    def setUp(self):
        """Set up test fixtures."""
        self.year = date.today().year
        self.cwd = tempfile.mkdtemp(prefix="TEST_", dir="/tmp")
        
        target = Female("FAM1", "F0", "001", "002", "003", target="1", age="30",
                       yob=str(self.year-30), cancers=Cancers(),
                       gtests=BWSGeneticTests.default_factory())
        self.pedigree = BwaPedigree(people=[target])
        self.pedigree.add_parents(target)

    def tearDown(self):
        """Clean up test fixtures."""
        TestCase.tearDown(self)
        shutil.rmtree(self.cwd)

    def test_write_batch_file_with_empty_calc_ages(self):
        """Test write_batch_file with empty calc_ages list."""
        ped_file = os.path.join(self.cwd, "test.ped")
        batch_file = os.path.join(self.cwd, "test_empty.bat")
        
        # Empty list should default to [0]
        result = self.pedigree.write_batch_file(
            pedigree_file_name=ped_file,
            filepath=batch_file,
            calc_ages=[],
            model_settings=settings.BC_MODEL
        )
        
        self.assertTrue(os.path.exists(result))

    def test_write_batch_file_age_range(self):
        """Test write_batch_file generates appropriate age range."""
        ped_file = os.path.join(self.cwd, "test.ped")
        batch_file = os.path.join(self.cwd, "test_range.bat")
        
        # Target is age 30, should calculate ages in appropriate range
        result = self.pedigree.write_batch_file(
            pedigree_file_name=ped_file,
            filepath=batch_file,
            model_settings=settings.BC_MODEL
        )
        
        self.assertTrue(os.path.exists(result))
        with open(result, 'r') as f:
            content = f.read()
            # Check that ages are calculated (should have multiple "9" entries for age calculations)
            self.assertGreater(content.count('9'), 0)

    def test_write_pedigree_file_target_attributes(self):
        """Test writing input pedigree file for fortran with target-specific attributes."""
        filepath = os.path.join(self.cwd, "test_target.ped")
        target = self.pedigree.get_target()
        target.target = "1"  # Ensure it's marked as target
        
        result = self.pedigree.write_pedigree_file(
            filepath=filepath,
            risk_factor_code="12345678",
            hgt=170.5,
            model_settings=settings.BC_MODEL
        )
        
        self.assertTrue(os.path.exists(result))

        with open(result, 'r') as f:
            content = f.read()
            self.assertIn("12345678", content)
            self.assertIn("170.5000", content)

    def test_write_pedigree_file_non_target_person(self):
        """Test write_pedigree_file with non-target person in pedigree."""
        filepath = os.path.join(self.cwd, "test_nontarget.ped")
        
        # Add a non-target person
        nontarget = Female("FAM1", "F1", "004", "0", "0", target="0", 
                          age="35", yob=str(self.year-35), cancers=Cancers())
        self.pedigree.people.append(nontarget)
        
        result = self.pedigree.write_pedigree_file(filepath=filepath,
                                                   model_settings=settings.BC_MODEL)
        
        self.assertTrue(os.path.exists(result))


class PedigreeGetMethods(TestCase):
    """Test pedigree get methods."""

    def setUp(self):
        """Set up test fixtures."""
        self.year = date.today().year
        
        target = Female("FAM1", "F0", "001", "002", "003", target="1", age="30",
                       yob=str(self.year-30), cancers=Cancers())
        self.pedigree = BwaPedigree(people=[target])
        self.pedigree.add_parents(target)

    def test_get_person_found(self):
        """Test getting a person by ID."""
        person = self.pedigree.get_person("001")
        self.assertIsNotNone(person)
        self.assertEqual(person.pid, "001")

    def test_get_person_by_name_found(self):
        """Test getting a person by name."""
        target = self.pedigree.get_target()
        target.name = "ProbandsName"
        person = self.pedigree.get_person_by_name("ProbandsName")
        self.assertIsNotNone(person)
        self.assertEqual(person.name, "ProbandsName")

    def test_get_siblings(self):
        """Test getting siblings."""
        target = self.pedigree.get_target()
        father = self.pedigree.get_person(target.fathid)
        mother = self.pedigree.get_person(target.mothid)
        
        # Add a sibling
        sibling = Female("FAM1", "F1", "004", father.pid, mother.pid, 
                        age="28", yob=str(self.year-28), cancers=Cancers())
        self.pedigree.people.append(sibling)
        
        siblings, siblings_yob = self.pedigree.get_siblings(target)
        self.assertGreater(len(siblings), 0)

    def test_get_twins(self):
        """Test getting twins."""
        target = self.pedigree.get_target()
        target.mztwin = "A"
        
        # Add a twin
        twin = Female("FAM1", "F1", "004", target.fathid, target.mothid,
                     age=target.age, yob=target.yob, mztwin="A",
                     cancers=Cancers())
        self.pedigree.people.append(twin)
        
        twins_dict = self.pedigree.get_twins()
        self.assertIn("A", twins_dict)
        self.assertEqual(len(twins_dict["A"]), 2)

    def test_unconnected_all_connected(self):
        """Test unconnected returns empty list when all are connected."""
        unconnected = self.pedigree.unconnected()
        self.assertEqual(len(unconnected), 0)

    def test_unconnected_has_unconnected(self):
        """Test unconnected detects unconnected members."""
        # Add an unconnected person
        unconnected_person = Female("FAM1", "F1", "100", "0", "0",
                                   age="50", yob="1960", cancers=Cancers())
        self.pedigree.people.append(unconnected_person)
        
        unconnected_ids = self.pedigree.unconnected()
        self.assertIn("100", unconnected_ids)
