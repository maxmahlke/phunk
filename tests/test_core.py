from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import phunk

# Observations from Gehrels 1956
PHASE = [0.57, 1.09, 3.20, 10.99, 14.69, 20.42]
MAG = [6.555, 6.646, 6.793, 7.130, 7.210, 7.414]

OBS_8988 = pd.read_csv(Path(__file__).parent / "data/Hansenkoharcheck.csv")


@pytest.mark.parametrize(
    "model,hmag_expected",
    [
        ("HG", 6.504386696308138),
        ("HG1G2", 6.435448791494396),
        ("HG12", 6.47026229130713),
        ("HG12S", 6.4186901600356965),
        ("sHG1G2", 6.435448793849064),
        ("LinExp", -3730.422348715857),
    ],
)
def test_fitting_models(model, hmag_expected):
    """Test fitting all models."""
    # User has mag and phase and fits all models
    phase = PHASE
    mag = MAG

    # For sHG1G2
    ra = [184.568391, 184.338311, 183.433022, 180.361071, 179.244605, 178.693809]
    dec = [-2.467547, -2.361879, -1.944443, -0.500051, 0.048319, 0.410682]

    pc = phunk.PhaseCurve(phase=phase, mag=mag, ra=ra, dec=dec)
    pc.fit(models=[model])

    np.testing.assert_almost_equal(
        getattr(getattr(pc, model), "H" if model != "LinExp" else "a"), hmag_expected
    )


def test_multi_band():
    """Test instantiation and fitting of multi-band phase curves."""

    # Using local cache of ATLAS phase curve of (8988) 1979 MA4
    pc = phunk.PhaseCurve(
        phase=OBS_8988.phase,
        mag=OBS_8988.mred,
        band=OBS_8988.band,
        ra=np.degrees(OBS_8988.ra),
        dec=np.degrees(OBS_8988.dec),
    )

    assert pc.bands == ["c", "o"]

    pc.fit(models=["HG1G2", "sHG1G2"])

    # ------
    # HG1G2

    # Assert that we have band specific parameters
    for param in ["Hc", "Ho", "G1c", "G1o", "G2c", "G2c"]:
        assert hasattr(pc.HG1G2, param)

    # Assert that we do not have general parameters
    for param in ["H", "G1", "G2"]:
        assert not hasattr(pc.HG1G2, param)

    # ------
    # sHG1G2

    # Assert that we have band specific parameters
    for param in ["Hc", "Ho", "G1c", "G1o", "G2c", "G2c"]:
        assert hasattr(pc.sHG1G2, param)

    # Assert that we do not have general parameters
    for param in ["H", "G1", "G2"]:
        assert not hasattr(pc.sHG1G2, param)

    # Except for the shape
    for param in ["R", "alpha", "delta"]:
        assert hasattr(pc.sHG1G2, param)
