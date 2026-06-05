"""
© 2023 University of Cambridge
SPDX-FileCopyrightText: 2023 University of Cambridge
SPDX-License-Identifier: GPL-3.0-or-later
"""

import importlib.util
import os
import sys
import tempfile
from pathlib import Path
from unittest import TestCase
from unittest.mock import mock_open, patch

import bws.settings as settings
from vcf2prs.exception import Vcf2PrsError


class SettingsTests(TestCase):

    def test_get_prs_alpha_dict_handles_string_and_dict_values(self):
        ''' Test that get_prs_alpha_dict correctly processes both string and dict formats for PRS alpha values. '''
        model = {
            'PRS_REFERENCE_FILES': {
                'EUROPEAN': [
                    ('BCAC 313', 'BCAC_313_PRS.prs'),
                    ('BCAC 307', {'alpha': 2.5})
                ],
                'AFRICAN': [],
                'EAST_ASIAN': [],
                'SOUTH_ASIAN': []
            }
        }

        alpha_dict = settings.get_prs_alpha_dict(model)
        self.assertEqual(alpha_dict, {'BCAC 313': 0.501, 'BCAC 307': 2.5})

    def test_get_prs_alpha_dict_rejects_duplicate_keys(self):
        ''' Test that get_prs_alpha_dict raises an error when duplicate keys are present in the PRS reference files. '''
        model = {
            'PRS_REFERENCE_FILES': {
                'EUROPEAN': [
                    ('BCAC 313', 'BCAC_313_PRS.prs'),
                    ('BCAC 313', 'BCAC_307_PRS.prs')
                ],
                'AFRICAN': [],
                'EAST_ASIAN': [],
                'SOUTH_ASIAN': []
            }
        }

        with self.assertRaises(AssertionError):
            settings.get_prs_alpha_dict(model)

    @patch('bws.settings.line_that_contain', 
           side_effect=UnicodeDecodeError('utf-8', b'', 0, 1, 'invalid'))
    @patch('builtins.open', new_callable=mock_open)
    def test_unicode_error_raises_vcf2prserror(self, mock_file, mock_line):
        """Test UnicodeDecodeError raises VCF2PrsError."""

        with self.assertRaisesRegex(Vcf2PrsError, "Unable to open the file"):
            settings.get_alpha("test.txt")

    def test_missing_prs_file_returns_empty_string(self):
        """Test that get_alpha returns an empty string when the PRS reference file is missing."""
        self.assertTrue(settings.get_alpha("test.txt") == "")

