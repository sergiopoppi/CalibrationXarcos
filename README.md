# Calibration with XARCOS

Calibrates XARCOS spectra stored in a CLASS file (from `discos2class`) to Jy.

For every spectrum, three factors are applied:

1. **Opacity**: `exp(tau / sin(el))`, which gives Ta*
2. **Elevation gain curve**: `1 / gain_curve(el)`, with `el` in degrees
3. **Gain**: `1 / gain(t)`, which gives Jy

`tau` and `gain` are **linearly interpolated** in time, per polarization (LCP/RCP),
from the input tables. Outside the tables' time range the edge value is used.
A warning is printed when the closest sample is more than `--max-dt` days away.

## Requirements

- GILDAS with `pyclass`
- Python 3, `numpy` (`pip install -r requirements.txt`)

## Usage

```bash
python src/calibrate_jy.py obs.class src/tau.txt src/gains.txt
```

Output: `obs_tstar_gain_jy.med`.

Options:

| Option | Description |
|---|---|
| `--keep-intermediate` | Runs three passes and also writes `obs_tstar.med` and `obs_tstar_gain.med` |
| `--max-dt DAYS` | Warning threshold for the distance to the closest tau/gain sample (default 0.25) |
| `--min-elevation DEG` | Spectra below this elevation are skipped (default 5) |
| `-v` | Logs the factor applied to every spectrum |

## Input file format

`tau.txt` and `gains.txt`: whitespace-separated columns.

```
MJD   LCP   RCP
```

The polarization of each spectrum comes from the CLASS line name, which must contain `LCP` or `RCP`.

## Tests

The computation logic is in `src/calibration_core.py`, which does not depend on pyclass:

```bash
pip install pytest
python -m pytest tests
```
