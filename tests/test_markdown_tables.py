"""
Tests for ``markdown_tables.py`` — MarkdownTable builder and shared helpers.

Covers MarkdownTable (rendering, alignment, escaping, edge cases),
best_in_group, fmt_val, is_failure, status_mark, eq_wins, compare_bases,
and the _read_log_info / _read_log_footer helpers.
"""

import os
import sys
import tempfile

import numpy as np
import pytest

CODE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
sys.path.insert(0, CODE_DIR)

from equivariant_pde_find.markdown_tables import (
    MarkdownTable,
    best_in_group,
    compare_bases,
    eq_wins,
    fmt_val,
    is_failure,
    status_mark,
    _read_log_footer,
    _read_log_info,
)


# ======================================================================
# MarkdownTable — construction
# ======================================================================


class TestMarkdownTableConstruction:
    """Construction and validation."""

    def test_default_alignment_is_left(self):
        t = MarkdownTable(["A", "B", "C"])
        assert t._align == ["left", "left", "left"]

    def test_custom_alignment(self):
        t = MarkdownTable(["A", "B"], align=["center", "right"])
        assert t._align == ["center", "right"]

    def test_alignment_length_mismatch_raises(self):
        with pytest.raises(ValueError, match="align has"):
            MarkdownTable(["A", "B"], align=["left"])

    def test_invalid_alignment_raises(self):
        with pytest.raises(ValueError, match="Invalid alignment"):
            MarkdownTable(["A"], align=["middle"])

    def test_add_row_wrong_length_raises(self):
        t = MarkdownTable(["A", "B"])
        with pytest.raises(ValueError, match="Expected 2 cells"):
            t.add_row("only_one")

    def test_add_row_accepts_non_string(self):
        t = MarkdownTable(["A"])
        t.add_row(42)
        t.add_row(3.14)


# ======================================================================
# MarkdownTable — rendering
# ======================================================================


class TestMarkdownTableRender:
    """Rendering correctness."""

    def test_empty_table(self):
        t = MarkdownTable(["Col1", "Col2"])
        rendered = t.render()
        lines = rendered.split("\n")
        assert len(lines) == 2  # header + separator
        assert lines[0] == "| Col1 | Col2 |"
        assert lines[1] == "|:-----|:-----|"

    def test_single_row(self):
        t = MarkdownTable(["Name", "Value"])
        t.add_row("foo", "bar")
        rendered = t.render()
        # "Value" (5 chars) is wider than "bar" (3 chars)
        assert "| foo  | bar   |" in rendered

    def test_left_alignment_separator(self):
        t = MarkdownTable(["A"], align=["left"])
        sep = t.render().split("\n")[1]
        assert sep == "|:--|"

    def test_center_alignment_separator(self):
        t = MarkdownTable(["A"], align=["center"])
        sep = t.render().split("\n")[1]
        assert sep == "|:-:|"

    def test_right_alignment_separator(self):
        t = MarkdownTable(["A"], align=["right"])
        sep = t.render().split("\n")[1]
        assert sep == "|--:|"

    def test_mixed_alignment(self):
        t = MarkdownTable(["L", "C", "R"], align=["left", "center", "right"])
        rendered = t.render()
        lines = rendered.split("\n")
        assert lines[0] == "| L | C | R |"
        assert lines[1] == "|:--|:-:|--:|"

    def test_column_widths_match_content(self):
        t = MarkdownTable(["A", "LongColumn"])
        t.add_row("x", "short")
        rendered = t.render()
        # "A" → width 1, "LongColumn" → width 10
        # Cell for "x": " x " (3 chars), cell for "short": " short      " (12 chars)
        assert "| x | short      |" in rendered

    def test_int_and_float_cells(self):
        t = MarkdownTable(["Int", "Float"])
        t.add_row(1, 3.14159)
        t.add_row(-5, 0.0)
        rendered = t.render()
        assert "| 1   | 3.14159 |" in rendered
        assert "| -5  | 0.0     |" in rendered

    def test_multiple_rows(self):
        t = MarkdownTable(["#", "Value"])
        for i in range(3):
            t.add_row(str(i), str(i * 10))
        rendered = t.render()
        lines = rendered.split("\n")
        assert len(lines) == 5  # header + sep + 3 rows

    def test_pipe_in_cell_is_not_escaped(self):
        """The MarkdownTable does NOT escape pipe characters in cells.
        This test documents the current behavior."""
        t = MarkdownTable(["Col"])
        t.add_row("a|b")
        rendered = t.render()
        # The pipe char appears literally — may break table rendering
        assert "a|b" in rendered


# ======================================================================
# fmt_val
# ======================================================================


class TestFmtVal:
    def test_float_default_precision(self):
        assert fmt_val(3.14159) == "3.1416"

    def test_float_custom_precision(self):
        assert fmt_val(1.0, d=2) == "1.00"

    def test_int_treated_as_float(self):
        assert fmt_val(5) == "5.0000"

    def test_string_passthrough(self):
        assert fmt_val("hello") == "hello"

    def test_none_passthrough(self):
        assert fmt_val(None) == "None"

    def test_nan_string(self):
        # "nan" string cannot be converted to float
        assert fmt_val("nan") == "nan"


# ======================================================================
# is_failure
# ======================================================================


class TestIsFailure:
    def test_exact_recovery_is_not_failure(self):
        r = {"exact_recovery": "1", "fp": "0", "fn": "0"}
        assert not is_failure(r)

    def test_fp_positive_is_failure(self):
        r = {"exact_recovery": "0", "fp": "2", "fn": "0"}
        assert is_failure(r)

    def test_fn_positive_is_failure(self):
        r = {"exact_recovery": "0", "fp": "0", "fn": "1"}
        assert is_failure(r)

    def test_both_positive_is_failure(self):
        r = {"exact_recovery": "0", "fp": "1", "fn": "3"}
        assert is_failure(r)

    def test_none_is_failure(self):
        assert is_failure(None)

    def test_missing_keys_is_failure(self):
        assert is_failure({})

    def test_nan_fp_is_failure(self):
        r = {"fp": "nan", "fn": "0"}
        assert is_failure(r)


# ======================================================================
# status_mark
# ======================================================================


class TestStatusMark:
    def test_exact(self):
        assert "\u2713" in status_mark({"exact_recovery": "1"})

    def test_non_exact(self):
        assert "\u2717" in status_mark({"exact_recovery": "0"})

    def test_none(self):
        assert status_mark(None) == ""


# ======================================================================
# best_in_group
# ======================================================================


class TestBestInGroup:
    def test_exact_beats_nonexact(self):
        rows = [
            {"exact_recovery": "0", "fp": "0", "fn": "0", "E_coeff": "0.001"},
            {"exact_recovery": "1", "fp": "0", "fn": "0", "E_coeff": "0.01"},
        ]
        best = best_in_group(rows)
        assert best["exact_recovery"] == "1"

    def test_lower_fp_beats_higher(self):
        rows = [
            {"exact_recovery": "0", "fp": "5", "fn": "0", "E_coeff": "0.001"},
            {"exact_recovery": "0", "fp": "2", "fn": "0", "E_coeff": "0.1"},
        ]
        best = best_in_group(rows)
        assert best["fp"] == "2"

    def test_lower_fn_breaks_tie(self):
        rows = [
            {"exact_recovery": "0", "fp": "2", "fn": "3", "E_coeff": "0.001"},
            {"exact_recovery": "0", "fp": "2", "fn": "1", "E_coeff": "0.1"},
        ]
        best = best_in_group(rows)
        assert best["fn"] == "1"

    def test_lower_E_breaks_fn_tie(self):
        rows = [
            {"exact_recovery": "0", "fp": "2", "fn": "1", "E_coeff": "0.05"},
            {"exact_recovery": "0", "fp": "2", "fn": "1", "E_coeff": "0.01"},
        ]
        best = best_in_group(rows)
        assert best["E_coeff"] == "0.01"

    def test_single_row(self):
        row = {"exact_recovery": "0", "fp": "3", "fn": "1", "E_coeff": "0.05"}
        assert best_in_group([row]) is row

    def test_nan_values_treated_as_inf(self):
        """NaN fp/fn are treated as infinity → lose to any finite value."""
        rows = [
            {"exact_recovery": "0", "fp": "nan", "fn": "nan", "E_coeff": "nan"},
            {"exact_recovery": "0", "fp": "1", "fn": "0", "E_coeff": "0.01"},
        ]
        best = best_in_group(rows)
        assert best["fp"] == "1"

    def test_missing_keys_default(self):
        """Missing exact_recovery defaults to '0' (not exact)."""
        rows = [
            {"fp": "0", "fn": "0", "E_coeff": "0.01"},
            {"exact_recovery": "1", "fp": "0", "fn": "0", "E_coeff": "0.1"},
        ]
        best = best_in_group(rows)
        assert best["exact_recovery"] == "1"


# ======================================================================
# eq_wins
# ======================================================================


class TestEqWins:
    def test_both_none(self):
        assert eq_wins(None, None) == ""

    def test_eq_exact_ne_fail(self):
        """eq exact (fp=0, fn=0), ne has fn>0 → eq clear win."""
        eq = {"exact_recovery": "1", "fp": "0", "fn": "0", "E_coeff": "0.01"}
        ne = {"exact_recovery": "0", "fp": "0", "fn": "1", "E_coeff": "0.1"}
        assert eq_wins(eq, ne) == "\U0001f7e2"  # eq clear win

    def test_ne_exact_eq_fail(self):
        """ne exact, eq has fp>0 → non-eq clear win."""
        eq = {"exact_recovery": "0", "fp": "1", "fn": "0", "E_coeff": "0.1"}
        ne = {"exact_recovery": "1", "fp": "0", "fn": "0", "E_coeff": "0.01"}
        assert eq_wins(eq, ne) == "\U0001f534"  # non-eq clear win

    def test_both_failed(self):
        """Both have fp>0 → both failures → em dash."""
        eq = {"exact_recovery": "0", "fp": "1", "fn": "0", "E_coeff": "0.1"}
        ne = {"exact_recovery": "0", "fp": "0", "fn": "1", "E_coeff": "0.2"}
        assert eq_wins(eq, ne) == "\u2014"  # both failed

    def test_both_exact_tie(self):
        """Both exact with same fp+fn and E within 1e-4 → tie."""
        eq = {"exact_recovery": "1", "fp": "0", "fn": "0", "E_coeff": "0.01"}
        ne = {"exact_recovery": "1", "fp": "0", "fn": "0", "E_coeff": "0.01"}
        assert eq_wins(eq, ne) == "\u26aa"  # tie

    def test_both_exact_eq_lower_E(self):
        """Both exact, same fp+fn=0, eq has lower E → eq marginal win."""
        eq = {"exact_recovery": "1", "fp": "0", "fn": "0", "E_coeff": "0.001"}
        ne = {"exact_recovery": "1", "fp": "0", "fn": "0", "E_coeff": "0.01"}
        assert eq_wins(eq, ne) == "\U0001f7e9"  # eq marginal

    def test_both_exact_ne_lower_E(self):
        """Both exact, same fp+fn=0, ne has lower E → non-eq marginal."""
        eq = {"exact_recovery": "1", "fp": "0", "fn": "0", "E_coeff": "0.01"}
        ne = {"exact_recovery": "1", "fp": "0", "fn": "0", "E_coeff": "0.001"}
        assert eq_wins(eq, ne) == "\U0001f7e5"  # non-eq marginal

    def test_eq_lower_fp_fn_both_not_exact(self):
        """Both not exact → both failures → em dash regardless of totals."""
        eq = {"exact_recovery": "0", "fp": "1", "fn": "0", "E_coeff": "0.1"}
        ne = {"exact_recovery": "0", "fp": "2", "fn": "1", "E_coeff": "0.01"}
        assert eq_wins(eq, ne) == "\u2014"


# ======================================================================
# compare_bases
# ======================================================================


class TestCompareBases:
    def test_modified_wins_exact(self):
        o = {"exact_recovery": "0", "fp": "0", "fn": "1", "E_coeff": "0.1"}
        m = {"exact_recovery": "1", "fp": "0", "fn": "0", "E_coeff": "0.01"}
        assert compare_bases(o, m) == "\U0001f7e2"

    def test_original_wins_exact(self):
        o = {"exact_recovery": "1", "fp": "0", "fn": "0", "E_coeff": "0.01"}
        m = {"exact_recovery": "0", "fp": "0", "fn": "1", "E_coeff": "0.1"}
        assert compare_bases(o, m) == "\U0001f534"

    def test_modified_lower_total(self):
        o = {"exact_recovery": "0", "fp": "2", "fn": "1", "E_coeff": "0.01"}
        m = {"exact_recovery": "0", "fp": "1", "fn": "0", "E_coeff": "0.1"}
        assert compare_bases(o, m) == "\U0001f7e2"

    def test_tie_lower_E_modified(self):
        o = {"exact_recovery": "0", "fp": "1", "fn": "1", "E_coeff": "0.1"}
        m = {"exact_recovery": "0", "fp": "1", "fn": "1", "E_coeff": "0.01"}
        assert compare_bases(o, m) == "\U0001f7e9"

    def test_tie_lower_E_original(self):
        o = {"exact_recovery": "0", "fp": "1", "fn": "1", "E_coeff": "0.01"}
        m = {"exact_recovery": "0", "fp": "1", "fn": "1", "E_coeff": "0.1"}
        assert compare_bases(o, m) == "\U0001f7e5"

    def test_tie_same_E(self):
        o = {"exact_recovery": "0", "fp": "1", "fn": "1", "E_coeff": "0.01"}
        m = {"exact_recovery": "0", "fp": "1", "fn": "1", "E_coeff": "0.01"}
        assert compare_bases(o, m) == "\u26aa"


# ======================================================================
# _read_log_info / _read_log_footer
# ======================================================================


class TestReadLogInfo:
    def test_no_log_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            csv_path = os.path.join(tmp, "nonexistent.csv")
            assert _read_log_info(csv_path) == {}
            assert _read_log_footer(csv_path) == ""

    def test_valid_log_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            csv_path = os.path.join(tmp, "sweep.csv")
            log_path = os.path.join(tmp, "sweep.log")
            with open(log_path, "w") as f:
                f.write("# comment line\n")
                f.write("total_wall_time_min: 12.5\n")
                f.write("omp_num_threads: 16\n")
                f.write("cpu_model: AMD EPYC\n")
            info = _read_log_info(csv_path)
            assert info["total_wall_time_min"] == "12.5"
            assert info["omp_num_threads"] == "16"
            assert info["cpu_model"] == "AMD EPYC"

    def test_log_footer_format(self):
        with tempfile.TemporaryDirectory() as tmp:
            csv_path = os.path.join(tmp, "sweep.csv")
            log_path = os.path.join(tmp, "sweep.log")
            with open(log_path, "w") as f:
                f.write("total_wall_time_min: 5.0\n")
                f.write("omp_num_threads: 8\n")
                f.write("cpu_model: Test CPU\n")
            footer = _read_log_footer(csv_path)
            assert "5.0 min" in footer
            assert "Test CPU" in footer
            assert "OMP_NUM_THREADS=8" in footer

    def test_empty_log(self):
        with tempfile.TemporaryDirectory() as tmp:
            csv_path = os.path.join(tmp, "sweep.csv")
            log_path = os.path.join(tmp, "sweep.log")
            with open(log_path, "w") as f:
                f.write("# only comments\n")
            assert _read_log_info(csv_path) == {}
