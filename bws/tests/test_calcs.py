"""
BOADICEA web-service testing.

© 2026 University of Cambridge
SPDX-FileCopyrightText: 2026 University of Cambridge
SPDX-License-Identifier: GPL-3.0-or-later
"""
import os

import pytest
from collections import OrderedDict
from unittest.mock import MagicMock, patch, mock_open
from django.test import TestCase, RequestFactory
from rest_framework.request import Request
from rest_framework.exceptions import ValidationError

from bws.exceptions import TimeOutException, ModelError
from bws.calc.calcs import Predictions
from bws.calc.model import ModelParams, ModelOpts
from bws.pedigree import Pedigree


# ---------------------------------------------------------------------------
# Helpers / fixtures
# ---------------------------------------------------------------------------

def make_mock_pedigree(size=5, siblings=None):
    ''' Return a minimal Pedigree mock. '''
    pedi = MagicMock(spec=Pedigree)
    pedi.people = [MagicMock() for _ in range(size)]
    target = MagicMock()
    pedi.get_target.return_value = target
    sib_list = siblings if siblings is not None else []
    pedi.get_siblings.return_value = (sib_list, [])
    return pedi


def make_mock_request(user_id=42):
    ''' Return a minimal Request mock with a user ID. '''
    rf = RequestFactory()
    request = rf.get("/")
    request.user = MagicMock()
    request.user.id = user_id
    return Request(request)


BC_MODEL_SETTINGS = {
    "NAME": "BC",
    "HOME": "/fake/home",
    "EXE": "boadicea",
    "INCIDENCE": "/fake/incidence/",
    "CALCS": ["carrier_probs", "remaining_lifetime"],
    "GENES": ["BRCA1", "BRCA2", "PALB2", "CHEK2", "ATM"],
}

OC_MODEL_SETTINGS = dict(BC_MODEL_SETTINGS, NAME="OC")
PC_MODEL_SETTINGS = dict(BC_MODEL_SETTINGS, NAME="PC")


# Minimal Fortran-like output used in multiple tests
SAMPLE_OUTPUT = """\

## LIFETIME RISK [option: -rl]
FollowUp Age,Censor Age,Risk
20,80,0.1200146

## 10-YEAR RISK  (NHS PROTOCOL) [option: -rj]
FollowUp Age,Censor Age,Risk
25,35,0.0023429
26,36,0.0028008
27,37,0.0033194
28,38,0.0039014
29,39,0.0045496
30,40,0.0052707

## 10-YEAR RISK [option: -ry]
FollowUp Age,Censor Age,Risk
40,50,0.0171806

## PROBABILITIES [option: -p]
NO_PATHOGENIC_VARIANTS,BRCA1,BRCA2,PALB2,CHEK2,ATM,BARD1,RAD51C,RAD51D
0.9821809,0.0012713,0.0020312,0.0012758,0.0074259,0.0035743,0.0008529,0.0006940,0.0006937

## REMAINING LIFETIME RISK [option: -rr]
Censor Age,Risk
26,0.0000735
27,0.0001709
28,0.0002957
29,0.0004508
30,0.0006415
35,0.0023429
40,0.0059089
45,0.0125627
50,0.0229952
55,0.0359022
60,0.0500926
65,0.0664033
70,0.0842939
75,0.1019159
80,0.1198730

"""


class TestGetNiceness(TestCase):
    ''' Tests for the _get_niceness function, which determines how "nice" the 
        BOADICEA process should be based on the pedigree size. The niceness is 
        calculated as the number of siblings divided by 15, capped at 19. 
        This is to prevent overloading the server with large pedigrees, while 
        allowing some flexibility for smaller ones. '''
    @pytest.mark.req_WS_CORE_100
    def test_no_siblings_small_pedigree(self):
        ''' A small pedigree with no siblings should have niceness 0. '''
        pedi = make_mock_pedigree(size=5)
        self.assertEqual(Predictions._get_niceness(pedi), 0)

    @pytest.mark.req_WS_CORE_100
    def test_no_siblings_large_pedigree(self):
        ''' A large pedigree with no siblings should have niceness 19 (capped). '''
        pedi = make_mock_pedigree(size=300)
        self.assertEqual(Predictions._get_niceness(pedi), 19)

    @pytest.mark.req_WS_CORE_100
    def test_no_siblings_medium_pedigree(self):
        ''' A medium pedigree with no siblings should have niceness 2. '''
        pedi = make_mock_pedigree(size=30)
        # 30 / 15 == 2
        self.assertEqual(Predictions._get_niceness(pedi), 2)

    @pytest.mark.req_WS_CORE_100
    def test_few_siblings(self):
        ''' A pedigree with a few siblings should have niceness equal to the number of siblings. '''
        siblings = [MagicMock() for _ in range(3)]
        pedi = make_mock_pedigree(siblings=siblings)
        self.assertEqual(Predictions._get_niceness(pedi), 3)

    @pytest.mark.req_WS_CORE_100
    def test_many_siblings_capped_at_19(self):
        ''' A pedigree with many siblings should have niceness capped at 19. '''
        siblings = [MagicMock() for _ in range(25)]
        pedi = make_mock_pedigree(siblings=siblings)
        self.assertEqual(Predictions._get_niceness(pedi), 19)

    @pytest.mark.req_WS_CORE_100
    def test_exactly_19_siblings(self):
        ''' A pedigree with exactly 19 siblings should have niceness 19. '''
        siblings = [MagicMock() for _ in range(19)]
        pedi = make_mock_pedigree(siblings=siblings)
        self.assertEqual(Predictions._get_niceness(pedi), 19)

    @pytest.mark.req_WS_CORE_100
    def test_custom_factor(self):
        ''' The niceness can be calculated using a custom factor instead of 15. '''
        pedi = make_mock_pedigree(size=20)
        # 20 / 10 == 2
        self.assertEqual(Predictions._get_niceness(pedi, factor=10), 2)


class TestGetVersion(TestCase):
    ''' Tests for the _get_version function, which calls the BOADICEA executable
        with a version flag and parses the output to determine the version string. 
        This is important for logging and debugging purposes, as different versions
        of BOADICEA may have different features or bugs. The tests mock the subprocess
        call to simulate different outputs and error conditions. '''

    @pytest.mark.req_WS_CORE_101
    @patch("bws.calc.calcs.Popen")
    def test_returns_version_string(self, mock_popen):
        ''' The function should return the version string output by the executable. '''
        proc = MagicMock()
        proc.communicate.return_value = (b"boadicea.exe model 7.3.2, version 0.6.0\n", b"")
        proc.wait.return_value = 0
        mock_popen.return_value = proc

        version = Predictions._get_version(model=BC_MODEL_SETTINGS)
        self.assertEqual(version, "boadicea model 7.3.2, version 0.6.0")

    @pytest.mark.req_WS_CORE_101
    @patch("bws.calc.calcs.Popen")
    def test_strips_exe_and_newline(self, mock_popen):
        ''' The function should strip the .exe suffix and any trailing newline from
            the version string. '''
        proc = MagicMock()
        proc.communicate.return_value = (b"bboadicea.exe model 7.3.2, version 0.6.0\n", b"")
        proc.wait.return_value = 0
        mock_popen.return_value = proc

        version = Predictions._get_version(model=BC_MODEL_SETTINGS)
        self.assertNotIn(".exe", version)
        self.assertNotIn("\n", version)

    @pytest.mark.req_WS_CORE_101
    @patch("bws.calc.calcs.Popen")
    def test_non_zero_exit_raises_model_error(self, mock_popen):
        ''' If the executable returns a non-zero exit code, the function should raise
            a ModelError with the error message. '''
        proc = MagicMock()
        proc.communicate.return_value = (b"", b"fatal error\n")
        proc.wait.return_value = 1
        mock_popen.return_value = proc

        with self.assertRaisesRegex(ModelError, "fatal error"):
            Predictions._get_version(model=BC_MODEL_SETTINGS)

    @pytest.mark.req_WS_CORE_101
    @patch("bws.calc.calcs.Popen")
    def test_timeout_raises_timeout_exception(self, mock_popen):
        ''' If the subprocess call times out, the function should raise a
            TimeOutException with an appropriate message. '''
        from subprocess import TimeoutExpired
        proc = MagicMock()
        proc.communicate.side_effect = TimeoutExpired(cmd="boadicea", timeout=30)
        mock_popen.return_value = proc

        with self.assertRaisesRegex(TimeOutException, "Request has timed out."):
            Predictions._get_version(model=BC_MODEL_SETTINGS)

    @pytest.mark.req_WS_CORE_101
    @patch("bws.calc.calcs.Popen")
    def test_unexpected_exception_propagates(self, mock_popen):
        ''' If the subprocess call raises an unexpected exception (e.g. 
            FileNotFoundError), it should propagate up and not be caught as
            a ModelError. '''
        mock_popen.side_effect = FileNotFoundError("exe not found")

        with self.assertRaisesRegex(FileNotFoundError, "exe not found"):
            Predictions._get_version(model=BC_MODEL_SETTINGS)


class TestRun(TestCase):
    ''' Tests for the run function, which is responsible for executing the BOADICEA risk calculation by calling the executable with the appropriate arguments and handling
        the output. The tests mock the subprocess call and file I/O to simulate
        different scenarios, such as a successful run, a run that returns an error,
        a run that times out, and a run with unsupported ethnicity for the prostate
        cancer model. The tests check that the function returns the expected output
        or raises the expected exceptions in each case. '''


    def _make_predictions_stub(self):
        """Return a bare Predictions instance without calling __init__."""
        return object.__new__(Predictions)

    def _base_model_params(self):
        mp = MagicMock(spec=ModelParams)
        mp.cancer_rates = "UK"
        mp.isashk = False
        mp.mutation_sensitivity = 1.0
        ethnicity = MagicMock()
        ethnicity.get_filename.return_value = "UK-pop.nml"
        mp.ethnicity = ethnicity
        return mp

    def _base_model_opts(self, out="can_risk.out"):
        opts = MagicMock(spec=ModelOpts)
        opts.out = out
        opts.get_cmd_line_opts.return_value = ["-rl", "-rr"]
        return opts

    @pytest.mark.req_WS_CORE_102
    @patch("bws.calc.calcs.open", new_callable=mock_open, read_data="output data")
    @patch("bws.calc.calcs.Popen")
    def test_successful_run_returns_file_contents(
        self, mock_popen, mock_file
    ):
        ''' When the subprocess call succeeds and produces output, the function
            should return the contents of the output file. '''
        proc = MagicMock()
        proc.communicate.return_value = (b"", b"")
        proc.wait.return_value = 0
        mock_popen.return_value = proc

        result = Predictions.run(
            request=make_mock_request(),
            bat_file="/tmp/risk.bat",
            model_opts=self._base_model_opts(),
            model_params=self._base_model_params(),
            cwd="/tmp",
            model=BC_MODEL_SETTINGS,
        )
        self.assertEqual(result, "output data")

    @pytest.mark.req_WS_CORE_102
    @patch("bws.calc.calcs.os.remove")
    @patch("bws.calc.calcs.open", new_callable=mock_open, read_data="output data")
    @patch("bws.calc.calcs.Popen")
    def test_successful_run_ignores_missing_previous_output_file(
        self, mock_popen, mock_file, mock_remove
    ):
        """If the previous output file does not exist, the run function should
        continue without raising and still return the new output."""
        proc = MagicMock()
        proc.communicate.return_value = (b"", b"")
        proc.wait.return_value = 0
        mock_popen.return_value = proc

        # os.remove raises OSError (e.g. file not found) — should be silently swallowed
        mock_remove.side_effect = OSError("no such file")

        model_opts = self._base_model_opts()   # out = e.g. "BC_predictions.txt"
        model_params = self._base_model_params()

        result = Predictions.run(
            request=make_mock_request(),
            bat_file="/tmp/risk.bat",
            model_opts=model_opts,
            model_params=model_params,
            cwd="/tmp",
            model=BC_MODEL_SETTINGS,
        )

        # The function should return the output despite the OSError on remove
        self.assertEqual(result, "output data")

        # Assert remove was called with the correct path — cwd + model_opts.out
        expected_out_path = os.path.join("/tmp", model_opts.out)
        mock_remove.assert_called_once_with(expected_out_path)

    @pytest.mark.req_WS_CORE_102
    @patch("bws.calc.calcs.os.remove")
    @patch("bws.calc.calcs.open", new_callable=mock_open, read_data="output data")
    @patch("bws.calc.calcs.Popen")
    def test_run_builds_correct_subprocess_command(
        self, mock_popen, mock_file, mock_remove
    ):
        """The run function should include param file and incidence file in the
        subprocess command line."""
        proc = MagicMock()
        proc.communicate.return_value = (b"", b"")
        proc.wait.return_value = 0
        mock_popen.return_value = proc

        model_params = self._base_model_params()
        opts = self._base_model_opts()

        result = Predictions.run(
            request=make_mock_request(),
            bat_file="/tmp/risk.bat",
            model_opts=opts,
            model_params=model_params,
            param_file="/tmp/params.params",
            cwd="/tmp",
            model=BC_MODEL_SETTINGS,
        )

        self.assertEqual(result, "output data")

        cmd = mock_popen.call_args[0][0]

        # Executable
        self.assertEqual(cmd[0], os.path.join(BC_MODEL_SETTINGS['HOME'], BC_MODEL_SETTINGS['EXE']))

        # Param file flag
        self.assertEqual(cmd[1], "-s")
        self.assertEqual(cmd[2], "/tmp/params.params")

        # Ethnicity flag and file
        self.assertEqual(cmd[3], "-e")
        expected_ethnicity = os.path.join(BC_MODEL_SETTINGS["HOME"], 'Data', 
                                          "coeffs-"+BC_MODEL_SETTINGS['NAME']+"_"+model_params.ethnicity.get_filename())
        self.assertEqual(cmd[4], expected_ethnicity)

        # bat file and incidence file
        self.assertEqual(cmd[-2], "/tmp/risk.bat")
        expected_incidence = BC_MODEL_SETTINGS['INCIDENCE'] + model_params.cancer_rates + ".nml"
        self.assertEqual(cmd[-1], expected_incidence)

    @pytest.mark.req_WS_CORE_102
    @patch("bws.calc.calcs.open", new_callable=mock_open, read_data="")
    @patch("bws.calc.calcs.Popen")
    def test_nonzero_exit_raises_model_error(
        self, mock_popen, mock_file
    ):
        ''' When the subprocess call returns a non-zero exit code, the function
            should raise a ModelError with the error message from the subprocess. '''
        proc = MagicMock()
        proc.communicate.return_value = (b"stdout", b"some fortran error\n")
        proc.wait.return_value = 1
        mock_popen.return_value = proc

        with self.assertRaisesRegex(ModelError, "some fortran error"):
            Predictions.run(
                request=make_mock_request(),
                bat_file="/tmp/risk.bat",
                model_opts=self._base_model_opts(),
                model_params=self._base_model_params(),
                cwd="/tmp",
                model=BC_MODEL_SETTINGS,
            )

    @pytest.mark.req_WS_CORE_102
    @patch("bws.calc.calcs.Popen")
    def test_timeout_raises_timeout_exception(
        self, mock_popen
    ):
        ''' When the subprocess call times out, the function should raise a
            TimeOutException with an appropriate message. '''
        from subprocess import TimeoutExpired
        proc = MagicMock()
        proc.communicate.side_effect = TimeoutExpired(cmd="boadicea", timeout=30)
        mock_popen.return_value = proc

        with self.assertRaisesRegex(TimeOutException, "Request has timed out."):
            Predictions.run(
                request=make_mock_request(),
                bat_file="/tmp/risk.bat",
                model_opts=self._base_model_opts(),
                model_params=self._base_model_params(),
                cwd="/tmp",
                model=BC_MODEL_SETTINGS,
            )

    @pytest.mark.req_WS_CORE_102
    @patch("bws.calc.calcs.open", new_callable=mock_open, read_data="")
    @patch("bws.calc.calcs.Popen")
    def test_pc_model_non_uk_ethnicity_raises(
        self, mock_popen, mock_file
    ):
        ''' When running the prostate cancer model, if the ethnicity specified in
            the model parameters is not UK, the function should raise a ModelError
            indicating that the prostate cancer model does not support UK ethnicity.
            This is because the prostate cancer model is currently only calibrated
            for UK populations, and using a different ethnicity would produce invalid
            results. The test mocks the model parameters to return a non-UK ethnicity
            and checks that the appropriate exception is raised. '''   
        mp = self._base_model_params()
        mp.ethnicity.get_filename.return_value = "UK-southAsian.nml"

        with self.assertRaisesRegex(ModelError, "prostate cancer model does not support UK ethnicity"):
            Predictions.run(
                request=make_mock_request(),
                bat_file="/tmp/risk.bat",
                model_opts=self._base_model_opts(),
                model_params=mp,
                cwd="/tmp",
                model=PC_MODEL_SETTINGS,
            )


class TestParseRisksOutput(TestCase):
    ''' Tests for the _parse_risks_output function, which takes the raw output
        from the BOADICEA executable and parses it into structured data for the
        different types of risk calculations (remaining lifetime risk, 10-year
        risk, NHS protocol risk) and mutation probabilities. The function should
        correctly identify the sections of the output corresponding to each type
        of calculation, parse the CSV data within those sections, and return it
        in a consistent format. The tests check that the parsing is done correctly
        for each type of risk and that disabled options return None. '''

    def _make_predictions(self, model_settings=None):
        p = object.__new__(Predictions)
        p.model_settings = model_settings or BC_MODEL_SETTINGS
        return p

    def _all_opts(self):
        opts = MagicMock(spec=ModelOpts)
        opts.rr = True
        opts.rl = True
        opts.ry = True
        opts.rj = True
        opts.probs = True
        return opts

    @pytest.mark.req_WS_CORE_103
    def test_remaining_lifetime_risk_parsed(self):
        ''' The remaining lifetime risk section of the output should be parsed into
        a list of dictionaries with age and risk values. '''
        p = self._make_predictions()
        _rl, rr, _ry, _rj, _mp = p._parse_risks_output(SAMPLE_OUTPUT, self._all_opts())
        self.assertEqual(len(rr), 15)
        self.assertEqual(rr[0]["age"], 26)
        self.assertEqual(rr[0]["breast cancer risk"]["decimal"], 0.0000735)
        self.assertEqual(rr[0]["breast cancer risk"]["percent"], 0)

    @pytest.mark.req_WS_CORE_103
    def test_lifetime_risk_parsed(self):
        ''' The lifetime risk section of the output should be parsed into a list of
        dictionaries with age and risk values. '''
        p = self._make_predictions()
        rl, _rr, _ry, _rj, _mp = p._parse_risks_output(SAMPLE_OUTPUT, self._all_opts())
        self.assertEqual(len(rl), 1)
        self.assertEqual(rl[0]["age"], 80)
        self.assertEqual(rl[0]["breast cancer risk"]["decimal"], 0.1200146)

    @pytest.mark.req_WS_CORE_103
    def test_10yr_risk_parsed(self):
        ''' The 10-year risk section of the output should be parsed into a list of
        dictionaries with age and risk values. '''
        p = self._make_predictions()
        _rl, _rr, ry, _rj, _mp = p._parse_risks_output(SAMPLE_OUTPUT, self._all_opts())
        self.assertEqual(len(ry), 1)
        self.assertEqual(ry[0]["age"], 50)
        self.assertEqual(ry[0]["breast cancer risk"]["decimal"], 0.0171806)

    @pytest.mark.req_WS_CORE_103
    def test_nhs_protocol_risk_parsed(self):
        ''' The NHS protocol risk section of the output should be parsed into a
        list of dictionaries with age and risk values. '''
        p = self._make_predictions()
        _rl, _rr, _ry, rj, _mp = p._parse_risks_output(SAMPLE_OUTPUT, self._all_opts())
        self.assertEqual(len(rj), 6)
        self.assertEqual(rj[0]["breast cancer risk"]["decimal"], 0.0023429)

    @pytest.mark.req_WS_CORE_103
    def test_mutation_probs_parsed(self):
        ''' The mutation probabilities section of the output should be parsed into a
        list of dictionaries with gene names and probability values. '''
        p = self._make_predictions()
        _rl, _rr, _ry, _rj, mp = p._parse_risks_output(SAMPLE_OUTPUT, self._all_opts())
        self.assertIsNotNone(mp)
        self.assertTrue(len(mp) > 0)
        self.assertIn("no mutation", mp[0])

    @pytest.mark.req_WS_CORE_103
    def test_disabled_opts_return_none(self):
        ''' Disabled risk options should return None. '''
        p = self._make_predictions()
        opts = MagicMock(spec=ModelOpts)
        opts.rr = False
        opts.rl = False
        opts.ry = False
        opts.rj = False
        opts.probs = False
        rl, rr, ry, rj, mp = p._parse_risks_output("", opts)
        self.assertIsNone(rl)
        self.assertIsNone(rr)
        self.assertIsNone(ry)
        self.assertIsNone(rj)
        self.assertIsNone(mp)

    @pytest.mark.req_WS_CORE_103
    def test_cancer_type_label_oc(self):
        ''' The cancer type label for ovarian cancer should be correctly
        identified. '''
        p = self._make_predictions(model_settings=OC_MODEL_SETTINGS)
        opts = MagicMock(spec=ModelOpts)
        opts.rr = True
        opts.rl = opts.ry = opts.rj = opts.probs = False
        output = "## REMAINING LIFETIME RISK\n,50,0.10\n"
        _rl, rr, _ry, _rj, _mp = p._parse_risks_output(output, opts)
        self.assertIn("ovarian cancer risk", rr[0])

    @pytest.mark.req_WS_CORE_103
    def test_cancer_type_label_pc(self):
        ''' The cancer type label for prostate cancer should be correctly
        identified. '''
        p = self._make_predictions(model_settings=PC_MODEL_SETTINGS)
        opts = MagicMock(spec=ModelOpts)
        opts.rr = True
        opts.rl = opts.ry = opts.rj = opts.probs = False
        output = "## REMAINING LIFETIME RISK\n,50,0.08\n"
        _rl, rr, _ry, _rj, _mp = p._parse_risks_output(output, opts)
        self.assertIn("prostate cancer risk", rr[0])

    @pytest.mark.req_WS_CORE_103
    def test_empty_output_returns_empty_arrays(self):
        ''' Empty output should return empty arrays. '''
        p = self._make_predictions()
        rl, rr, ry, rj, mp = p._parse_risks_output("", self._all_opts())
        self.assertEqual(rl, [])
        self.assertEqual(rr, [])
        self.assertEqual(ry, [])
        self.assertEqual(rj, [])
        self.assertEqual(mp, [])

    @pytest.mark.req_WS_CORE_103
    def test_parse_risks_output_ignores_comments_and_blank_lines(self):
        ''' The parser should ignore comment lines and blank lines within the
            output and still return structured risk arrays. '''
        p = self._make_predictions()
        output = (
            "# comment line\n"
            "## LIFETIME RISK\n"
            "FollowUp Age,Censor Age,Risk\n"
            "20,80,0.1200146\n"
            "\n"
            "## REMAINING LIFETIME RISK\n"
            "Censor Age,Risk\n"
            "26,0.0000735\n"
        )
        rl, rr, ry, rj, mp = p._parse_risks_output(output, self._all_opts())
        self.assertEqual(len(rl), 1)
        self.assertEqual(len(rr), 1)
        self.assertEqual(ry, [])
        self.assertEqual(rj, [])
        self.assertEqual(mp, [])


class TestParseProbsOutput(TestCase):
    ''' Tests for the _parse_probs_output function, which takes the raw output
    from the BOADICEA executable and parses the mutation probabilities section
    into a structured format. The function should correctly identify the
    probabilities section, parse the CSV data, and return a list of dictionaries
    with gene names and probability values in both decimal and percent formats.
    The tests check that the parsing is done correctly, that the "no mutation"
    probability is handled appropriately, that gene columns match the model
    settings, and that empty or malformed input is handled gracefully. '''

    BC_PROBS = "NO_PATHOGENIC_VARIANTS,BRCA1,BRCA2,PALB2,CHEK2,ATM,BARD1,RAD51C,RAD51D\n0.9821809,0.0012713,0.0020312,0.0012758,0.0074259,0.0035743,0.0008529,0.0006940,0.0006937"

    @pytest.mark.req_WS_CORE_104
    def test_no_mutation_probability_first(self):
        result = Predictions._parse_probs_output(self.BC_PROBS, BC_MODEL_SETTINGS)
        assert "no mutation" in result[0]
        self.assertEqual(result[0]["no mutation"]["decimal"], 0.9821809)
        self.assertEqual(result[0]["no mutation"]["percent"], 98.22)

    @pytest.mark.req_WS_CORE_104
    def test_gene_probabilities_present(self):
        result = Predictions._parse_probs_output(self.BC_PROBS, BC_MODEL_SETTINGS)
        genes = [list(r.keys())[0] for r in result[1:]]
        self.assertEqual(genes, BC_MODEL_SETTINGS["GENES"])

    @pytest.mark.req_WS_CORE_104
    def test_decimal_and_percent_fields(self):
        result = Predictions._parse_probs_output(self.BC_PROBS, BC_MODEL_SETTINGS)
        for entry in result:
            val = list(entry.values())[0]
            self.assertTrue("decimal" in val)
            self.assertTrue("percent" in val)
            self.assertEqual(val["percent"], round(val["decimal"] * 100, 2))

    @pytest.mark.req_WS_CORE_104
    def test_empty_probs_returns_empty_list(self):
        result = Predictions._parse_probs_output("", BC_MODEL_SETTINGS)
        self.assertEqual(result, [])

    @pytest.mark.req_WS_CORE_104
    def test_gene_column_mismatch_raises(self):
        bad_probs = "BRCA1,WRONG_GENE,PALB2,CHEK2,ATM\n0.90,0.02,0.03,0.02,0.02,0.01\n"
        with self.assertRaisesRegex(AssertionError, "RESULTS COLUMN MISMATCH"):
            Predictions._parse_probs_output(bad_probs, BC_MODEL_SETTINGS)

    @pytest.mark.req_WS_CORE_104
    def test_pc_model_uses_dynamic_gene_columns(self):
        pc_probs = "NO_PATHOGENIC_VARIANTS,BRCA1,BRCA2,HOXB13\n0.90,0.05,0.025,0.025\n"
        result = Predictions._parse_probs_output(pc_probs, PC_MODEL_SETTINGS)
        self.assertTrue("no mutation" in result[0])
        genes = [list(r.keys())[0] for r in result[1:]]
        self.assertEqual(genes, ["BRCA1", "BRCA2", "HOXB13"])


class TestIsCalculate(TestCase):
    ''' Tests for the is_calculate function, which checks if a given calculation is
        included in the list of calculations to be performed, e.g. "carrier_probs",
        "remaining_lifetime", "lifetime", "ten_year". '''

    @pytest.mark.req_WS_CORE_105
    def _make(self, calcs):
        p = object.__new__(Predictions)
        p.calcs = calcs
        return p

    @pytest.mark.req_WS_CORE_105
    def test_empty_calcs_always_true(self):
        p = self._make([])
        self.assertTrue(p.is_calculate("carrier_probs"))
        self.assertTrue(p.is_calculate("remaining_lifetime"))
        self.assertTrue(p.is_calculate("lifetime"))
        self.assertTrue(p.is_calculate("ten_year"))

    @pytest.mark.req_WS_CORE_105
    def test_calc_in_list_returns_true(self):
        p = self._make(["carrier_probs"])
        self.assertTrue(p.is_calculate("carrier_probs"))

    @pytest.mark.req_WS_CORE_105
    def test_calc_not_in_list_returns_false(self):
        p = self._make(["carrier_probs"])
        self.assertFalse(p.is_calculate("remaining_lifetime"))


class TestPredictionsInit(TestCase):
    ''' Tests for the __init__ method of the Predictions class, which initializes a
    Predictions instance with a pedigree, optional calculations to perform, and
    other parameters. The tests check that the method correctly validates the
    input parameters, sets default values for calculations based on model settings,
    and handles edge cases such as invalid calculation names or incorrect types for
    pedigree and model parameters. The tests use mocking to isolate the
    initialization logic from the actual risk calculation logic. '''

    @pytest.mark.req_WS_CORE_106
    def _patch_run_risks(self):
        return patch.object(Predictions, "_run_risks", return_value=None)

    def _patch_version_and_niceness(self):
        return (
            patch.object(Predictions, "_get_version", return_value="5.0"),
            patch.object(Predictions, "_get_niceness", return_value=0),
        )

    @pytest.mark.req_WS_CORE_106
    @patch("bws.calc.calcs.settings")
    def test_invalid_calc_raises_validation_error(self, mock_settings):
        ''' Test mocks the settings to define allowed calculations and checks that
        providing an invalid name results in the expected exception. '''
        mock_settings.BC_MODEL = BC_MODEL_SETTINGS
        mock_settings.ALLOWED_CALCS = ["carrier_probs", "remaining_lifetime"]
        pedi = make_mock_pedigree()
        with self.assertRaisesRegex(ValidationError, "Unknown calculation requested: invalid_calc"):
            Predictions(pedi, calcs=["invalid_calc"], run_risks=False,
                        model_settings=BC_MODEL_SETTINGS)

    @pytest.mark.req_WS_CORE_106
    @patch("bws.calc.calcs.settings")
    def test_pedigree_type_assertion(self, mock_settings):
        ''' Test mocks the settings to define allowed calculations and checks that
        providing a non-Pedigree object results in the expected exception. '''
        mock_settings.BC_MODEL = BC_MODEL_SETTINGS
        mock_settings.ALLOWED_CALCS = ["carrier_probs", "remaining_lifetime"]
        with self.assertRaisesRegex(AssertionError, "'not_a_pedigree' is not a Pedigree"):
            Predictions("not_a_pedigree", run_risks=False,
                        model_settings=BC_MODEL_SETTINGS)

    @pytest.mark.req_WS_CORE_106
    @patch("bws.calc.calcs.settings")
    def test_model_params_type_assertion(self, mock_settings):
        ''' Test mocks the settings to define allowed calculations and checks that
        providing a non-ModelParams object results in the expected exception. '''
        mock_settings.BC_MODEL = BC_MODEL_SETTINGS
        mock_settings.ALLOWED_CALCS = ["carrier_probs", "remaining_lifetime"]
        pedi = make_mock_pedigree()
        with self.assertRaisesRegex(AssertionError, "'not_model_params' is not a ModelParams"):
            Predictions(pedi, model_params="not_model_params", run_risks=False,
                        model_settings=BC_MODEL_SETTINGS)

    @pytest.mark.req_WS_CORE_106
    @patch("bws.calc.calcs.settings")
    def test_risk_factor_code_coerced_to_str(self, mock_settings):
        ''' Test mocks the settings to define allowed calculations and checks that
        providing a non-string risk factor code results in it being coerced to a
        string. '''
        mock_settings.BC_MODEL = BC_MODEL_SETTINGS
        mock_settings.ALLOWED_CALCS = ["carrier_probs", "remaining_lifetime"]
        pedi = make_mock_pedigree()
        with self._patch_run_risks():
            pred = Predictions(pedi, risk_factor_code=7, run_risks=False,
                               model_settings=BC_MODEL_SETTINGS)
        self.assertEqual(pred.risk_factor_code, "7")

    @pytest.mark.req_WS_CORE_106
    @patch("bws.calc.calcs.settings")
    def test_calcs_defaults_to_model_settings(self, mock_settings):
        ''' Test mocks the settings to define allowed calculations and checks that
        if no calculations are provided, it defaults to the model settings. '''
        mock_settings.BC_MODEL = BC_MODEL_SETTINGS
        mock_settings.ALLOWED_CALCS = BC_MODEL_SETTINGS["CALCS"]
        pedi = make_mock_pedigree()
        with self._patch_run_risks():
            pred = Predictions(pedi, run_risks=False, model_settings=BC_MODEL_SETTINGS)
        self.assertEqual(pred.calcs, BC_MODEL_SETTINGS["CALCS"])


class TestRunRisksIntegration(TestCase):
    """
    Integration-style tests for _run_risks using mocked _run_risk results.
    These test the branching logic in _run_risks without hitting the filesystem.
    """

    @pytest.mark.req_WS_CORE_107
    def _base_pred(self):
        p = object.__new__(Predictions)
        p.pedi = make_mock_pedigree()
        p.model_settings = BC_MODEL_SETTINGS
        p.model_params = MagicMock(spec=ModelParams)
        p.request = make_mock_request()
        p.cwd = "/tmp"
        p.risk_factor_code = "0"
        p.hgt = -1
        p.mdensity = None
        p.prs = None
        p.calcs = BC_MODEL_SETTINGS["CALCS"]
        return p

    @pytest.mark.req_WS_CORE_107
    @patch.object(Predictions, "_get_version", return_value="5.0")
    @patch.object(Predictions, "_get_niceness", return_value=0)
    @patch.object(Predictions, "_run_risk")
    def test_lifetime_risk_set_on_pred(self, mock_run_risk, _niceness, _version):
        ''' When the risk output includes lifetime risk (rl), it should be set as
        an attribute on the Predictions instance. '''
        p = self._base_pred()
        rl = [OrderedDict([("age", 80), ("breast cancer risk", {"decimal": 0.32, "percent": 32.0})])]
        # Primary risk call returns (rl, rr, ry, rj, mp)
        mock_run_risk.return_value = (rl, None, None, None, None)
        p._run_risks()
        self.assertTrue(hasattr(p, "lifetime_cancer_risk"))
        self.assertEqual(p.lifetime_cancer_risk, rl)    

    @pytest.mark.req_WS_CORE_107
    @patch.object(Predictions, "_get_version", return_value="5.0")
    @patch.object(Predictions, "_get_niceness", return_value=0)
    @patch.object(Predictions, "_run_risk")
    def test_mutation_probabilities_set(self, mock_run_risk, _niceness, _version):
        ''' When the risk output includes mutation probabilities (mp), they should
        be set as an attribute on the Predictions instance. '''
        p = self._base_pred()
        mp = [{"no mutation": {"decimal": 0.9, "percent": 90.0}}]
        mock_run_risk.return_value = (None, None, None, None, mp)
        p._run_risks()
        self.assertTrue(hasattr(p, "mutation_probabilties"))
        self.assertEqual(p.mutation_probabilties, mp)

    @pytest.mark.req_WS_CORE_107
    @patch.object(Predictions, "_get_version", return_value="5.0")
    @patch.object(Predictions, "_get_niceness", return_value=0)
    @patch.object(Predictions, "_run_risk")
    def test_baseline_risks_set_when_rr_returned(self, mock_run_risk, _niceness, _version):
        ''' When the risk output includes remaining lifetime risk (rr), a second
        call to _run_risk should be made to get the baseline risks, and they should
        be set as an attribute on the Predictions instance. '''
        p = self._base_pred()
        rr = [OrderedDict([("age", 55), ("breast cancer risk", {"decimal": 0.10, "percent": 10.0})])]
        baseline_rr = [OrderedDict([("age", 55), ("breast cancer risk", {"decimal": 0.08, "percent": 8.0})])]

        # First call (Risk): returns rr only
        # Second call (RemainingLifetimeBaselineRisk): returns baseline rr
        # Third call (RiskBaseline): returns rl=None, ry=None
        mock_run_risk.side_effect = [
            (None, rr, None, None, None),
            (None, baseline_rr, None, None, None),
            (None, None, None, None, None),
        ]
        p._run_risks()
        self.assertTrue(hasattr(p, "baseline_cancer_risks"))
        self.assertEqual(p.baseline_cancer_risks, baseline_rr)

    @pytest.mark.req_WS_CORE_107
    @patch.object(Predictions, "_get_version", return_value="5.0")
    @patch.object(Predictions, "_get_niceness", return_value=0)
    @patch.object(Predictions, "_run_risk")
    def test_nhs_protocol_falls_back_to_ry(self, mock_run_risk, _niceness, _version):
        ''' When rj is empty but ry is available, ten_yr_nhs_protocol should equal
            ry. This tests the fallback logic in _run_risks for the NHS protocol risk.

            ry: 10-yr risk (proband's age set at 40y, censor age set at 50y)
            rj: 10-yr risk (NHS protocol for young women at high risk)
        '''
        p = self._base_pred()
        rl = [{"age": 80}]
        ry = [{"age": 50}]
        mock_run_risk.side_effect = [
            (rl, None, ry, [], None),   # rj=[] triggers fallback
            (None, None, None, None, None),
            (None, None, None, None, None),
        ]
        p._run_risks()
        self.assertTrue(hasattr(p, "ten_yr_nhs_protocol"))
        self.assertEqual(p.ten_yr_nhs_protocol, ry)

    @pytest.mark.req_WS_CORE_107
    @patch.object(Predictions, "_get_version", return_value="5.0")
    @patch.object(Predictions, "_get_niceness", return_value=0)
    @patch.object(Predictions, "_run_risk")
    def test_version_stored(self, mock_run_risk, _niceness, _version):
        ''' The version string returned by _get_version should be stored as an
        attribute on the Predictions instance. '''
        p = self._base_pred()
        mock_run_risk.return_value = (None, None, None, None, None)
        p._run_risks()
        self.assertEqual(p.version, "5.0")
