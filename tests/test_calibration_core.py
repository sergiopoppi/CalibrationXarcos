import math
import os
import sys

import numpy as np
import pytest

SRC = os.path.join(os.path.dirname(__file__), "..", "src")
sys.path.insert(0, SRC)

import calibration_core as cc  # noqa: E402


def test_read_table_sorts_and_shapes(tmp_path):
    f = tmp_path / "t.txt"
    f.write_text("2.0 0.2 0.3\n1.0 0.1 0.15\n")
    t = cc.read_table(str(f))
    assert t.shape == (2, 3)
    assert list(t[:, 0]) == [1.0, 2.0]


def test_read_table_rejects_bad_files(tmp_path):
    f = tmp_path / "t.txt"
    f.write_text("1.0 0.1\n2.0 0.2\n")
    with pytest.raises(ValueError):
        cc.read_table(str(f))


def test_repo_data_files_load():
    assert cc.read_table(os.path.join(SRC, "tau.txt")).shape == (6, 3)
    assert cc.read_table(os.path.join(SRC, "gains.txt")).shape == (5, 3)


def test_scan_mjd_matches_original_formula():
    dobs, ut = -2541.0, math.radians(90.0)  # 6h UT
    # original: Time(dobs + 2460549.5 + ut_fraction, format='jd').mjd
    expected = dobs + 2460549.5 + 0.25 - 2400000.5
    assert cc.scan_mjd(dobs, ut) == pytest.approx(expected)


@pytest.mark.parametrize("line,pol", [("H2O-RCP", "RCP"), ("h2o_lcp", "LCP")])
def test_polarization_from_line(line, pol):
    assert cc.polarization_from_line(line) == pol


@pytest.mark.parametrize("line", ["H2O", "RCP-LCP"])
def test_polarization_from_line_rejects_ambiguous(line):
    with pytest.raises(ValueError):
        cc.polarization_from_line(line)


TABLE = np.array([[10.0, 0.1, 1.0], [11.0, 0.3, 2.0]])


def test_interpolate_midpoint_per_polarization():
    assert cc.interpolate(TABLE, 10.5, "LCP") == pytest.approx(0.2)
    assert cc.interpolate(TABLE, 10.5, "RCP") == pytest.approx(1.5)


def test_interpolate_clamps_and_warns_when_far(caplog):
    assert cc.interpolate(TABLE, 20.0, "LCP", max_dt=0.25) == pytest.approx(0.3)
    assert "days away" in caplog.text


def test_interpolate_no_warning_when_close(caplog):
    cc.interpolate(TABLE, 10.1, "LCP", max_dt=0.25)
    assert caplog.text == ""


def test_gain_curve_matches_original_quadratic():
    c = cc.GAIN_CURVE_POLY
    for el in (10.0, 45.0, 80.0):
        assert cc.gain_curve(el) == pytest.approx(c[0] + c[1] * el + c[2] * el ** 2)


def test_opacity_factor_zenith():
    assert cc.opacity_factor(0.1, math.pi / 2) == pytest.approx(math.exp(0.1))
