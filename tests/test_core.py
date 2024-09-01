import numpy as np

import phunk

# Observations from Gehrels 1956
PHASE = [0.57, 1.09, 3.20, 10.99, 14.69, 20.42]
MAG = [6.555, 6.646, 6.793, 7.130, 7.210, 7.414]


def test_fitting_models():
    """Test fitting all different models."""


def test_phase_target_ephemerides():
    """Test providing phase angle versus target and ephermerides."""

    # User has mag and phase and fits all models
    phase = PHASE
    mag = MAG

    pc = phunk.PhaseCurve(phase=phase, mag=mag)
    pc.fit(models=phunk.models.MODELS)

    assert np.testing.assert_almost_equal(pc.HG1G2.H, 6.430901601343539)

    # User has mag, phase, ra, dec
    # -> Can fit HG1G2, sHG1G2

    # user has mag, epoch, target
    # -> Can fit HG1G2, sHG1G2
