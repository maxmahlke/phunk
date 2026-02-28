import numpy as np


def cos_aspect_angle(ra, dec, ra0, dec0):
    """Compute the cosine of the aspect angle

    This angle is computed from the coordinates of the target and
    the coordinates of its pole.
    See Eq 12.4 "Introduction to Ephemerides and Astronomical Phenomena", IMCCE

    Parameters
    ----------
    ra: float
        Right ascension of the target in radians.
    dec: float
        Declination of the target in radians.
    ra0: float
        Right ascension of the pole in radians.
    dec0: float
        Declination of the pole in radians.

    Returns
    -------
    float: The cosine of the aspect angle
    """
    return np.sin(dec) * np.sin(dec0) + np.cos(dec) * np.cos(dec0) * np.cos(ra - ra0)


def rotation_phase(t, W0, W1, t0):
    """Compute the rotational phase

    This angle is computed from the location of the prime meridian at
    at reference epoch (W0, t0), and an angular velocity (W1)
    See Eq 12.1 "Introduction to Ephemerides and Astronomical Phenomena", IMCCE

    Parameters
    ----------
    t: float
        Time (JD)
    W0: float
        Location of the prime meridian at reference epoch (radian)
    W1: float
        Angular velocity of the target in radians/day.
    t0: float
        Reference epoch (JD)

    Returns
    -------
    float: The rotational phase W (radian)
    """
    return W0 + W1 * (t - t0)


def subobserver_longitude(ra, dec, ra0, dec0, W):
    """Compute the subobserver longitude (radian)

    This angle is computed from the coordinates of the target,
    the coordinates of its pole, and its rotation phase
    See Eq 12.4 "Introduction to Ephemerides and Astronomical Phenomena", IMCCE

    Parameters
    ----------
    ra: float
        Right ascension of the target in radians.
    dec: float
        Declination of the target in radians.
    ra0: float
        Right ascension of the pole in radians.
    dec0: float
        Declination of the pole in radians.
    W: float
        Rotation phase of the target in radians.

    Returns
    -------
    float: The subobserver longitude in radians.
    """
    x = -np.cos(dec0) * np.sin(dec) + np.sin(dec0) * np.cos(dec) * np.cos(ra - ra0)
    y = -(np.cos(dec) * np.sin(ra - ra0))
    return W - np.arctan2(x, y)
