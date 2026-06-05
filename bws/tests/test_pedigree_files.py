"""
Pedigree file writing tests.

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

from django.test import TestCase
from django.conf import settings

from bws.cancer import Cancers, Cancer, BWSGeneticTests, CanRiskGeneticTests, GeneticTest, PathologyTests, PathologyTest
from bws.pedigree import BwaPedigree, CanRiskPedigree
from bws.person import Male, Female
from bws.risk_factors.ethnicity import ONSEthnicity
from bws.risk_factors.mdensity import Volpara, Birads, Stratus
from collections import namedtuple

# Simple PRS object
PRS = namedtuple('PRS', ['alpha', 'zscore'])


class PedigreeFilesTests(TestCase):
    """Test pedigree file writing methods."""

    def setUp(self):
        """Set up test fixtures."""
        self.year = date.today().year
        self.cwd = tempfile.mkdtemp(prefix="TEST_", dir="/tmp")
        
        # Create basic pedigree
        target = Female("FAM1", "F0", "001", "002", "003", target="1", age="30",
                       yob=str(self.year-30), cancers=Cancers(),
                       gtests=BWSGeneticTests.default_factory())
        self.pedigree = BwaPedigree(people=[target])
        self.pedigree.add_parents(target)
        
        # Create CanRisk pedigree
        target_cr = Female("FAM2", "F0", "001", "002", "003", target="1", age="40",
                          yob=str(self.year-40), cancers=Cancers(),
                          gtests=CanRiskGeneticTests.default_factory())
        self.canrisk_pedigree = CanRiskPedigree(people=[target_cr])
        self.canrisk_pedigree.add_parents(target_cr, gtests=CanRiskGeneticTests.default_factory())

    def tearDown(self):
        """Clean up test fixtures."""
        TestCase.tearDown(self)
        shutil.rmtree(self.cwd)

    def test_write_pedigree_file_basic(self):
        """Test writing input pedigree file for fortran with default parameters."""
        filepath = os.path.join(self.cwd, "test.ped")
        result = self.pedigree.write_pedigree_file(filepath=filepath, 
                                                   model_settings=settings.BC_MODEL)
        
        self.assertTrue(os.path.exists(result))
        with open(result, 'r') as f:
            lines = f.readlines()
            self.assertGreater(len(lines), 0)
            # Check header line
            self.assertIn("I3", lines[0])

    def test_write_pedigree_file_with_risk_factor(self):
        """Test writing input pedigree file for fortran with risk factor code."""
        filepath = os.path.join(self.cwd, "test_rf.ped")
        result = self.pedigree.write_pedigree_file(
            filepath=filepath,
            risk_factor_code="12345678",
            model_settings=settings.BC_MODEL
        )
        
        self.assertTrue(os.path.exists(result))
        with open(result, 'r') as f:
            lines = f.readlines()
            self.assertGreater(len(lines), 2)

    def test_write_pedigree_file_oc_model(self):
        """Test writing input pedigree file for fortran with OC model."""
        filepath = os.path.join(self.cwd, "test_oc.ped")
        target = self.pedigree.get_target()
        target.cancers = Cancers(bc1=Cancer(), bc2=Cancer(), oc=Cancer("35"), 
                                prc=Cancer(), pac=Cancer())
        
        result = self.pedigree.write_pedigree_file(
            filepath=filepath,
            model_settings=settings.OC_MODEL
        )
        
        self.assertTrue(os.path.exists(result))

    def test_write_pedigree_file_pc_model(self):
        """Test writing input pedigree file for fortran with PC model."""
        filepath = os.path.join(self.cwd, "test_pc.ped")
        # Create a male prostate cancer patient
        target = Male("FAM3", "M0", "001", "002", "003", target="1", age="50",
                     yob=str(self.year-50), cancers=Cancers(),
                     gtests=BWSGeneticTests.default_factory())
        pedigree_pc = BwaPedigree(people=[target])
        pedigree_pc.add_parents(target)
        
        result = pedigree_pc.write_pedigree_file(
            filepath=filepath,
            model_settings=settings.PC_MODEL
        )
        
        self.assertTrue(os.path.exists(result))

    def test_write_pedigree_file_with_height(self):
        """Test writing input pedigree file for fortran with height parameter."""
        filepath = os.path.join(self.cwd, "test_hgt.ped")
        result = self.pedigree.write_pedigree_file(
            filepath=filepath,
            hgt=165.5,
            model_settings=settings.BC_MODEL
        )
        
        self.assertTrue(os.path.exists(result))

    def test_write_pedigree_file_with_birads(self):
        """Test writing input pedigree file for fortran with Birads mammographic density."""
        filepath = os.path.join(self.cwd, "test_birads.ped")
        birads = Birads("1")
        result = self.pedigree.write_pedigree_file(
            filepath=filepath,
            mdensity=birads,
            model_settings=settings.BC_MODEL
        )
        
        self.assertTrue(os.path.exists(result))

    def test_write_pedigree_file_with_volpara(self):
        """Test writing input pedigree file for fortran with Volpara mammographic density."""
        filepath = os.path.join(self.cwd, "test_volpara.ped")
        volpara = Volpara("15.0")
        result = self.pedigree.write_pedigree_file(
            filepath=filepath,
            mdensity=volpara,
            model_settings=settings.BC_MODEL
        )
        
        self.assertTrue(os.path.exists(result))

    def test_write_pedigree_file_with_stratus(self):
        """Test writing input pedigree file for fortran with Stratus mammographic density."""
        filepath = os.path.join(self.cwd, "test_stratus.ped")
        stratus = Stratus("0.8")
        result = self.pedigree.write_pedigree_file(
            filepath=filepath,
            mdensity=stratus,
            model_settings=settings.BC_MODEL
        )
        
        self.assertTrue(os.path.exists(result))

    def test_write_pedigree_file_with_prs(self):
        """Test writing input pedigree file for fortran with polygenic risk score."""
        filepath = os.path.join(self.cwd, "test_prs.ped")
        prs = PRS(alpha=0.15, zscore=1.5)
        result = self.pedigree.write_pedigree_file(
            filepath=filepath,
            prs=prs,
            model_settings=settings.BC_MODEL
        )
        
        self.assertTrue(os.path.exists(result))

    def test_write_pedigree_file_invalid_mdensity(self):
        """Test writing input pedigree file for fortran with invalid mammographic density type."""
        filepath = os.path.join(self.cwd, "test_invalid_md.ped")
        
        with self.assertRaises(Exception) as cm:
            self.pedigree.write_pedigree_file(
                filepath=filepath,
                mdensity="invalid",
                model_settings=settings.BC_MODEL
            )
        self.assertIn("Unsupported", str(cm.exception))

    def test_write_param_file(self):
        """Test writing model parameter file."""
        filepath = os.path.join(self.cwd, "params")
        result = self.pedigree.write_param_file(
            filepath=filepath,
            model_settings=settings.BC_MODEL
        )
        
        self.assertTrue(os.path.exists(result))
        with open(result, 'r') as f:
            content = f.read()
            self.assertIn("&settings", content)
            self.assertIn("POPULATION_ALLELE_FRQ", content)

    def test_write_param_file_ashkenazi(self):
        """Test writing model parameter file for Ashkenazi population."""
        filepath = os.path.join(self.cwd, "params_ash")
        result = self.pedigree.write_param_file(
            filepath=filepath,
            isashk=True,
            model_settings=settings.BC_MODEL
        )
        
        self.assertTrue(os.path.exists(result))
        with open(result, 'r') as f:
            content = f.read()
            self.assertIn("PEDIGREE_ALLELE_FRQ", content)

    def test_write_batch_file_default_ages(self):
        """Test writing batch file with default calculated ages."""
        ped_file = os.path.join(self.cwd, "test.ped")
        batch_file = os.path.join(self.cwd, "test.bat")
        
        result = self.pedigree.write_batch_file(
            pedigree_file_name=ped_file,
            filepath=batch_file,
            model_settings=settings.BC_MODEL
        )
        
        self.assertTrue(os.path.exists(result))
        with open(result, 'r') as f:
            lines = f.readlines()
            self.assertGreater(len(lines), 0)

    def test_write_batch_file_custom_ages(self):
        """Test writing batch file with custom ages."""
        ped_file = os.path.join(self.cwd, "test.ped")
        batch_file = os.path.join(self.cwd, "test_custom.bat")
        
        custom_ages = [35, 40, 45, 50]
        result = self.pedigree.write_batch_file(
            pedigree_file_name=ped_file,
            filepath=batch_file,
            calc_ages=custom_ages,
            model_settings=settings.BC_MODEL
        )
        
        self.assertTrue(os.path.exists(result))
        with open(result, 'r') as f:
            content = f.read()
            self.assertIn("9", content)

    def test_write_batch_file_single_age(self):
        """Test writing batch file with single integer age."""
        ped_file = os.path.join(self.cwd, "test.ped")
        batch_file = os.path.join(self.cwd, "test_single.bat")
        
        result = self.pedigree.write_batch_file(
            pedigree_file_name=ped_file,
            filepath=batch_file,
            calc_ages=45,
            model_settings=settings.BC_MODEL
        )
        
        self.assertTrue(os.path.exists(result))

    def test_write_batch_file_pc_model(self):
        """Test writing batch file with PC model."""
        # Create male target
        target = Male("FAM3", "M0", "001", "002", "003", target="1", age="50",
                     yob=str(self.year-50), cancers=Cancers(),
                     gtests=BWSGeneticTests.default_factory())
        pedigree_pc = BwaPedigree(people=[target])
        pedigree_pc.add_parents(target)
        
        ped_file = os.path.join(self.cwd, "test_pc.ped")
        batch_file = os.path.join(self.cwd, "test_pc.bat")
        
        result = pedigree_pc.write_batch_file(
            pedigree_file_name=ped_file,
            filepath=batch_file,
            model_settings=settings.PC_MODEL
        )
        
        self.assertTrue(os.path.exists(result))

    def test_get_columns(self):
        """Test getting column headers."""
        columns = self.pedigree.get_columns()
        self.assertIsNotNone(columns)
        self.assertGreater(len(columns), 0)
        self.assertEqual(columns[0], "FamID")

    def test_get_column_idx(self):
        """Test getting column index by name."""
        idx = BwaPedigree.get_column_idx("FamID")
        self.assertEqual(idx, 0)
        
        idx = BwaPedigree.get_column_idx("Name")
        self.assertEqual(idx, 1)

    def test_canrisk_get_column_idx_canrisk1(self):
        """Test getting column index for CanRisk v1."""
        idx = CanRiskPedigree.get_column_idx("FamID", file_type="canrisk1")
        self.assertEqual(idx, 0)

    def test_canrisk_get_column_idx_canrisk2(self):
        """Test getting column index for CanRisk v2."""
        idx = CanRiskPedigree.get_column_idx("BRIP1", file_type="canrisk2")
        self.assertEqual(idx, 25)
        idx = CanRiskPedigree.get_column_idx("ER:PR:HER2:CK14:CK56", file_type="canrisk2")
        self.assertEqual(idx, 26)

    def test_canrisk_get_column_idx_canrisk4(self):
        """Test getting column index for CanRisk v4."""
        idx = CanRiskPedigree.get_column_idx("HOXB13", file_type="canrisk4")
        self.assertEqual(idx, 26)
        idx = CanRiskPedigree.get_column_idx("ER:PR:HER2:CK14:CK56", file_type="canrisk4")
        self.assertEqual(idx, 27)

    def test_canrisk_get_prs_bc(self):
        """Test getting BC PRS from CanRisk pedigree."""
        prs = PRS(alpha=0.1, zscore=1.2)
        target = self.canrisk_pedigree.get_target()
        self.canrisk_pedigree.bc_prs = prs
        
        result = self.canrisk_pedigree.get_prs('BC')
        self.assertEqual(result, prs)

    def test_canrisk_get_prs_oc(self):
        """Test getting OC PRS from CanRisk pedigree."""
        prs = PRS(alpha=0.15, zscore=1.5)
        self.canrisk_pedigree.oc_prs = prs
        
        result = self.canrisk_pedigree.get_prs('OC')
        self.assertEqual(result, prs)

    def test_canrisk_get_prs_pc(self):
        """Test getting PC PRS from CanRisk pedigree."""
        prs = PRS(alpha=0.2, zscore=2.0)
        self.canrisk_pedigree.pc_prs = prs
        
        result = self.canrisk_pedigree.get_prs('PC')
        self.assertEqual(result, prs)

    def test_canrisk_get_prs_not_set(self):
        """Test getting PRS when not set returns None."""
        result = self.canrisk_pedigree.get_prs('BC')
        self.assertIsNone(result)

    def test_canrisk_get_rfcode_bc(self):
        """Test getting BC risk factor code from CanRisk pedigree."""
        self.canrisk_pedigree.bc_risk_factor_code = "12345678"
        
        result = self.canrisk_pedigree.get_rfcode('BC')
        self.assertEqual(result, "12345678")

    def test_canrisk_get_rfcode_oc(self):
        """Test getting OC risk factor code from CanRisk pedigree."""
        self.canrisk_pedigree.oc_risk_factor_code = "87654321"
        
        result = self.canrisk_pedigree.get_rfcode('OC')
        self.assertEqual(result, "87654321")

    def test_canrisk_get_rfcode_not_set(self):
        """Test getting risk factor code when not set returns 0."""
        result = self.canrisk_pedigree.get_rfcode('BC')
        self.assertEqual(result, 0)


class PedigreeInitializationTests(TestCase):
    """Test pedigree initialization edge cases."""

    def setUp(self):
        """Set up test fixtures."""
        self.year = date.today().year

    def test_pedigree_with_canrisk_parameters(self):
        """Test creating pedigree with CanRisk specific parameters."""
        target = Female("FAM1", "F0", "001", "002", "003", target="1", age="35",
                       yob=str(self.year-35), cancers=Cancers(),
                       gtests=CanRiskGeneticTests.default_factory())
        
        pedigree = CanRiskPedigree(
            people=[target],
            file_type='canrisk2',
            hgt=170,
            mdensity=Birads("2"),
            bc_risk_factor_code="12345678",
            oc_risk_factor_code="87654321",
            bc_prs=PRS(alpha=0.1, zscore=1.0),
            oc_prs=PRS(alpha=0.15, zscore=1.5)
        )
        
        self.assertEqual(pedigree.hgt, 170)
        self.assertIsNotNone(pedigree.mdensity)
        self.assertEqual(pedigree.mdensity.get_display_str(), "BI-RADS 2")
        self.assertEqual(pedigree.bc_risk_factor_code, "12345678")
        self.assertEqual(pedigree.oc_risk_factor_code, "87654321")

    def test_pedigree_with_ethnicity(self):
        """Test creating pedigree with ethnicity information."""
        # Note: Ethnicity requires specific implementations like ONSEthnicity
        # For now, test that the pedigree can be created with ethnicity parameters
        target = Female("FAM1", "F0", "001", "002", "003", target="1", age="35",
                       yob=str(self.year-35), cancers=Cancers(),
                       gtests=CanRiskGeneticTests.default_factory())
        
        pedigree = CanRiskPedigree(
            people=[target],
            file_type='canrisk2',
            ons_ethnicity=ONSEthnicity("White", "English/Welsh/Scottish/Northern Irish/British"),
        )
        
        # Verify pedigree was created successfully
        self.assertIsNotNone(pedigree)

    def test_add_parents_with_specified_ids(self):
        """Test adding parents when parent IDs are already specified."""
        target = Female("FAM1", "F0", "001", "100", "101", target="1", age="30",
                       yob=str(self.year-30), cancers=Cancers())
        pedigree = BwaPedigree(people=[target])
        
        # Add parents with specified IDs
        father, mother = pedigree.add_parents(target)
        
        self.assertEqual(father.pid, "100")
        self.assertEqual(mother.pid, "101")
        self.assertEqual(len(pedigree.people), 3)

    def test_add_parents_without_specified_ids(self):
        """Test adding parents when parent IDs are not specified."""
        target = Female("FAM1", "F0", "001", "0", "0", target="1", age="30",
                       yob=str(self.year-30), cancers=Cancers())
        pedigree = BwaPedigree(people=[target])
        
        # Add parents with auto-generated IDs
        father, mother = pedigree.add_parents(target)
        
        self.assertNotEqual(father.pid, "0")
        self.assertNotEqual(mother.pid, "0")
        self.assertTrue(father.pid.startswith("001"))
        self.assertTrue(mother.pid.startswith("001"))
        self.assertEqual(len(pedigree.people), 3)
