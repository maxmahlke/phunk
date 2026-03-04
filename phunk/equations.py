import numpy as np
from sbpy import photometry as phot

from phunk.geometry import cos_aspect_angle, rotation_phase, subobserver_longitude


def func_shg1g2(pha, h, g1, g2, R, alpha, delta):
    """Return f(H, G1, G2, R, alpha, delta) part of the lightcurve in mag space

    Parameters
    ----------
    pha: array-like [3, N]
        List containing [phase angle in radians, RA in radians, Dec in radians]
    h: float
        Absolute magnitude in mag
    G1: float
        G1 parameter (no unit)
    G2: float
        G2 parameter (no unit)
    R: float
        Oblateness (no units)
    alpha: float
        RA of the spin (radian)
    delta: float
        Dec of the spin (radian)

    Returns
    ----------
    out: array of floats
        H - 2.5 log(f(G1G2)) - 2.5 log(f(R, spin))
    """
    ph = pha[0]
    ra = pha[1]
    dec = pha[2]

    # Standard HG1G2 part: h + f(alpha, G1, G2)
    func1 = phot.HG1G2().evaluate(ph, h, g1, g2)

    # Spin part
    geo = cos_aspect_angle(ra, dec, alpha, delta)
    func2 = 1 - (1 - R) * np.abs(geo)
    func2 = 2.5 * np.log10(func2)

    return func1 + func2


def func_socca(pha, h, g1, g2, alpha0, delta0, period, a_b, a_c, phi0):
    """Return f(H, G1, G2, alpha0, delta0, period, a_b, a_c, phi0) part of the lightcurve in mag space

    Notes
    -----
    Absolute magnitude is computed according to Ostro & Connelly (1984)

    Parameters
    ----------
    pha: array-like [6, N]
        List containing [phase angle in radians, RA in radians, Dec in radians, time (jd), RA sun in radians, DEC sun in radians ]
    h: float
        Absolute magnitude in mag
    G1: float
        G1 parameter (no unit)
    G2: float
        G2 parameter (no unit)
    alpha0: float
        RA of the spin (radian)
    delta0: float
        Dec of the spin (radian)
    period: float
        Sidereal rotation period (days)
    a_b: float
        Equatorial axes ratio
    a_c: float
        Polar axes ratio
    phi0: float
        Initial rotation phase at reference time t0 (radian)
    t0: float
        Reference time (jd)

    Notes
    -----
    Input times must be corrected from the light travel time,
    that is jd_lt = jd - d_obs / c_speed

    Returns
    -------
    out: array of floats
        H - 2.5 log(f(G1G2)) - 2.5 log(f(spin, shape))
    """
    ph = pha[0]
    ra = pha[1]
    dec = pha[2]
    ep = pha[3]
    ra_s = pha[4]
    dec_s = pha[5]

    # TBD: For the time being, we fix the reference time
    # Time( '2022-01-01T00:00:00', format='isot', scale='utc').jd
    # Kinda middle of ZTF
    # TODO: take the middle jd?
    t0 = ep.mean()

    # Standard HG1G2 part: h + f(alpha, G1, G2)
    func1 = phot.HG1G2().evaluate(ph, h, g1, g2)

    # Rotation
    W = rotation_phase(ep, phi0, 2 * np.pi / period, t0)

    # Sub-Earth (e.TQe):
    # Spin part
    cos_aspect = cos_aspect_angle(ra, dec, alpha0, delta0)
    cos_aspect_2 = cos_aspect**2
    sin_aspect_2 = 1 - cos_aspect_2

    # Sidereal
    rot_phase = subobserver_longitude(ra, dec, alpha0, delta0, W)

    # https://www.sciencedirect.com/science/article/pii/0019103588901261
    eQe = (
        sin_aspect_2 * (np.cos(rot_phase) ** 2 + (a_b**2) * np.sin(rot_phase) ** 2)
        + cos_aspect_2 * a_c**2
    )

    # Sub-Solar (s.TQs):

    cos_aspect_s = cos_aspect_angle(ra_s, dec_s, alpha0, delta0)
    cos_aspect_s_2 = cos_aspect_s**2
    sin_aspect_s_2 = 1 - cos_aspect_s_2

    # Sidereal
    rot_phase_s = subobserver_longitude(ra_s, dec_s, alpha0, delta0, W)

    sQs = (
        sin_aspect_s_2
        * (np.cos(rot_phase_s) ** 2 + (a_b**2) * np.sin(rot_phase_s) ** 2)
        + cos_aspect_s_2 * a_c**2
    )

    # Cross-term (e.TQs):

    # sin(Lamda), sin(Lamda_sun):
    sin_aspect = np.sqrt(sin_aspect_2)
    sin_aspect_s = np.sqrt(sin_aspect_s_2)

    eQs = (
        sin_aspect * np.cos(rot_phase) * sin_aspect_s * np.cos(rot_phase_s)
        + sin_aspect * np.sin(rot_phase) * sin_aspect_s * np.sin(rot_phase_s) * (a_b**2)
        + cos_aspect * cos_aspect_s * a_c**2
    )

    # Full spin-shape term:
    I_tot = (np.sqrt(eQe) + eQs / np.sqrt(sQs)) / 2
    I_tot = -2.5 * np.log10(I_tot)

    return func1 + I_tot


def residual_shg1g2(pars, phase, mag, weights, bands, ra, dec):
    """Build the system of equations to solve using the HG1G2 + spin model

    ```
    x = [
        R, alpha, delta,
        h_g, g_1_g, g_2_g,
        h_r, g_1_r, g_2_r
    ]
    ```

    """
    R, alpha, delta = pars["R"], pars["alpha"], pars["delta"]
    # filternames = np.unique(bands)

    # params = x[3:]
    params = [value for name, value in pars.items() if not name.startswith("delta")]
    params += [delta]
    params = params[:-3]
    filternames = [
        param.name[1:] for param in pars.values() if param.name.startswith("H")
    ]
    # params = list(pars.values())[:-3]
    nparams = len(params) / len(filternames)
    # assert int(nparams) == nparams, "You need to input all parameters for all bands"

    params_per_band = np.reshape(params, (len(filternames), int(nparams)))
    eqs = []
    for index, filtername in enumerate(filternames):
        mask = bands == filtername

        myfunc = (
            func_shg1g2(
                np.vstack(
                    [phase[mask].tolist(), ra[mask].tolist(), dec[mask].tolist()]
                ),
                params_per_band[index][0],
                params_per_band[index][1],
                params_per_band[index][2],
                R,
                alpha,
                delta,
            )
            - mag[mask]
        ) / (1 / weights[mask])  # weighting by inverse mag_err

        eqs = np.concatenate((eqs, myfunc))

    return np.ravel(eqs)


def residual_socca(pars, phase, mag, weights, bands, ra, dec, ep, ra_s, dec_s):
    """Build the system of equations to solve using the HG1G2 + spin model

    ```
    x = [
        R, alpha, delta,
        h_g, g_1_g, g_2_g,
        h_r, g_1_r, g_2_r
    ]
    ```

    """
    alpha, delta, period, a_b, a_c, W0 = (
        pars["alpha"],
        pars["delta"],
        pars["period"],
        pars["a_b"],
        pars["a_c"],
        pars["W0"],
    )

    # filternames = np.unique(bands)

    # params = x[3:]
    params = [value for name, value in pars.items()]
    # params += [delta]
    params = params[:-3]
    filternames = [
        param.name[1:] for param in pars.values() if param.name.startswith("H")
    ]
    # params = list(pars.values())[:-3]
    nparams = len(params) / len(filternames)
    # assert int(nparams) == nparams, "You need to input all parameters for all bands"

    params_per_band = np.reshape(params, (len(filternames), int(nparams)))
    eqs = []
    for index, filtername in enumerate(filternames):
        mask = bands == filtername

        myfunc = (
            func_socca(
                np.vstack(
                    [
                        phase[mask].tolist(),
                        ra[mask].tolist(),
                        dec[mask].tolist(),
                        ep[mask].tolist(),
                        ra_s[mask].tolist(),
                        dec_s[mask].tolist(),
                    ]
                ),
                params_per_band[index][0],
                params_per_band[index][1],
                params_per_band[index][2],
                alpha,
                delta,
                period,
                a_b,
                a_c,
                W0,
            )
            - mag[mask]
        )

        eqs = np.concatenate((eqs, myfunc))

    return np.ravel(eqs)

