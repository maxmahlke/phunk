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
        # ("sHG1G2", 6.430901601343539),
        ("LinExp", -3730.422348715857),
    ],
)
def test_fitting_models(model, hmag_expected):
    """Test fitting all models."""
    # User has mag and phase and fits all models
    phase = PHASE
    mag = MAG

    pc = phunk.PhaseCurve(phase=phase, mag=mag)
    pc.fit(models=[model])

    np.testing.assert_almost_equal(
        getattr(getattr(pc, model), "H" if model != "LinExp" else "a"), hmag_expected
    )


def test_multi_band():
    """Test instantiation and fitting of multi-band phase curves."""
    pc = phunk.PhaseCurve(phase=OBS_8988.phase, mag=OBS_8988.mred, band=OBS_8988.band)

    assert pc.bands == ["c", "o"]


# different bands -> when fitting any model, parameters get band suffix
