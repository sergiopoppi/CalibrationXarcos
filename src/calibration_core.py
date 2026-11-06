"""Pure-numpy calibration helpers for XARCOS spectra.

Nothing in this module depends on GILDAS/pyclass, so it can be imported
and tested on any machine.
"""

import logging
import math

import numpy as np

log = logging.getLogger(__name__)

# CLASS stores the observation date (``dobs``) as an integer day count.
# The original script converted it with ``JD = dobs + 2460549.5 + ut``,
# i.e. ``MJD = dobs + 60549 + ut`` (JD - MJD = 2400000.5).
CLASS_DOBS_MJD_OFFSET = 60549

# Elevation gain curve, coefficients in ascending order of power:
# gain(el) = c0 + c1*el + c2*el**2 + c3*el**3, with el in degrees.
GAIN_CURVE_POLY = (0.7929185, 0.005900533, -4.203179e-5, 0.0)

# Column of each polarization in tau.txt and gains.txt (column 0 is MJD).
POLARIZATION_COLUMN = {"LCP": 1, "RCP": 2}

# Default limits.
DEFAULT_MAX_DT_DAYS = 0.25
DEFAULT_MIN_ELEVATION_DEG = 5.0


def read_table(path):
    """Read a whitespace-separated ``MJD LCP RCP`` table, sorted by MJD."""
    data = np.atleast_2d(np.genfromtxt(path, dtype=float))
    if data.ndim != 2 or data.shape[1] < 3:
        raise ValueError("%s: expected at least 3 columns (MJD LCP RCP)" % path)
    if np.isnan(data).any():
        raise ValueError("%s: contains missing or non-numeric values" % path)
    return data[np.argsort(data[:, 0])]


def scan_mjd(dobs, ut_rad):
    """MJD of a scan from CLASS ``dobs`` (days) and ``ut`` (radians)."""
    return dobs + CLASS_DOBS_MJD_OFFSET + ut_rad / (2.0 * math.pi)


def polarization_from_line(line):
    """Return ``'RCP'`` or ``'LCP'`` from a CLASS line name."""
    line = str(line).upper()
    has_rcp, has_lcp = "RCP" in line, "LCP" in line
    if has_rcp == has_lcp:
        raise ValueError("cannot determine polarization from line name %r" % line)
    return "RCP" if has_rcp else "LCP"


def interpolate(table, mjd, polarization, max_dt=DEFAULT_MAX_DT_DAYS, label="value"):
    """Linearly interpolate the table column for ``polarization`` at ``mjd``.

    Outside the table's time range the edge value is used. A warning is
    logged when the closest sample is more than ``max_dt`` days away.
    """
    times = table[:, 0]
    values = table[:, POLARIZATION_COLUMN[polarization]]
    nearest_dt = float(np.min(np.abs(times - mjd)))
    if nearest_dt > max_dt:
        log.warning(
            "%s %s at MJD %.5f: closest sample is %.3f days away (limit %.3f)",
            label, polarization, mjd, nearest_dt, max_dt,
        )
    return float(np.interp(mjd, times, values))


def gain_curve(el_deg, coeffs=GAIN_CURVE_POLY):
    """Relative gain at elevation ``el_deg`` (degrees)."""
    return float(np.polynomial.polynomial.polyval(el_deg, coeffs))


def opacity_factor(tau, el_rad):
    """Atmospheric correction ``exp(tau / sin(el))`` (plane-parallel airmass)."""
    return math.exp(tau / math.sin(el_rad))
