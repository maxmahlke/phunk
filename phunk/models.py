"""phunk: Implementation of photometric models"""

import lmfit
import numpy as np
from sbpy import photometry as phot

MODELS = ["Linexp", "HG", "HG12", "HG12S", "HG1G2", "sHG1G2"]


class HG:
    """HG model from Bowell+ 1989"""


class Linexp:
    """Linear-Exponential Modeling from Muinonen+ 2002"""

    def __init__(self, a=15.0, b=0.15, d=0.15, k=1):
        self.a = a
        self.b = b
        self.d = d
        self.k = k

    def evaluate(self, phase, a=None, b=None, d=None, k=None):
        """Evaluate photometric phase curve model."""
        phase = np.radians(phase)

        if a is None:
            a = self.a
        if b is None:
            b = self.b
        if d is None:
            d = self.d
        if k is None:
            k = self.k

        return a * np.exp(-phase / d) + b + k * phase


class HG12S:
    """HG12* model from Penttilä+ 2016."""


class HG1G2:
    """HG1G2 model from Muinonenü 2010."""

    def __init__(self, H=15.0, G1=0.15, G2=0.15):
        self.name = "HG1G2"
        self.H = H
        self.G1 = G1
        self.G2 = G2

    def evaluate(self, phase, H=None, G1=None, G2=None):
        """H,G1,G2 phase curve model."""
        phase = np.radians(phase)

        if H is None:
            H = self.H
        if G1 is None:
            G1 = self.G1
        if G2 is None:
            G2 = self.G2

        return phot.HG1G2().evaluate(phase, H, G1, G2)

    def is_fittable(self, pc):
        """Check if phase curve can be fit."""

        # We need at least magnitudes
        assert (
            pc.mag is not None
        ), f"No magnitudes provided. Cannot fit {self.name} model."

        # If we have matching phase angles, all is good
        if pc.phase is not None:
            assert (
                len(pc.mag) == len(pc.phase)
            ), f"Provided 'phase' and 'mag' have different lengths: {len(pc.phase)} vs {len(pc.mag)}"
            return True

        # If no phase: We need epochs and a target name
        assert (
            pc.epoch is not None
        ), "If no 'phase' is provided you have to define the 'epoch'."
        assert (
            pc.target is not None
        ), "If no 'phase' is provided you have to define the 'target'."

        pc.compute_ephemerides()
        return True

    def fit(self, phase, mag, mag_err):
        """Fit a phase curve using the H, G1, G2 model from Muinonen+ 2010.

        Parameters
        ----------
        phase : list or np.ndarray
            The phase angles of observation in degree.
        mag : list or np.ndarray
            The reduced magnitudes.
        """

        # The model function
        model = lmfit.Model(phot.HG1G2.evaluate)

        # The fit parameters
        params = lmfit.Parameters()
        params.add("h", value=15, min=0, max=30)
        params.add("g1", value=0.15, min=0, max=1.0)

        # Add delta to implement inequality constraint
        # https://lmfit.github.io/lmfit-py/constraints.html#using-inequality-constraints
        params.add("delta", min=0, value=0.5, max=1, vary=True)
        params.add("g2", expr="delta-g1", min=0.0, max=1)

        # Ensure that errors are well behaved
        # mag_err[mag_err == 0] = np.nanmean(mag_err)
        weights = weights_from_phase(phase)  # used for weighting

        # And fit
        result = model.fit(
            mag,
            params,
            ph=np.radians(phase),
            method="least_squares",
            fit_kws={
                "loss": "soft_l1",
            },
            weights=weights,
        )

        self.H = result.params["h"].value
        self.G1 = result.params["g1"].value
        self.G2 = result.params["g2"].value
        self.H_err = result.params["h"].stderr
        self.G1_err = result.params["g1"].stderr
        self.G2_err = result.params["g2"].stderr


class sHG1G2:
    """sHG1G2 phase curve model."""

    def __init__(self, obs, H=15.0, G1=0.15, G2=0.15):
        self.name = "sHG1G2"
        self.obs = obs
        self.H = H
        self.G1 = G1
        self.G2 = G2

    def eval(self, band):
        """sHG1G2 phase curve model. Adapted implementation from J. Peloton.

        Return f(H, G1, G2, R, alpha0, delta0) part of the lightcurve in mag space

        Parameters
        ----------
        pha: array-like [3, N]
            List containing [phase angle in degrees, RA in radians, Dec in radians]
        h: float
            Absolute magnitude in mag
        G1: float
            G1 parameter (no unit)
        G2: float
            G2 parameter (no unit)
        R: float
            Oblateness (no units)
        alpha0: float
            RA of the spin (radian)
        delta0: float
            Dec of the spin (radian)

        Returns
        ----------
        out: array of floats
            H - 2.5 log(f(G1G2)) - 2.5 log(f(R, spin))
        """
        obs_ = self.obs[self.obs.band == band]

        # Standard HG1G2 part: h + f(alpha, G1, G2)
        func1 = phot.HG1G2().evaluate(
            np.radians(obs_.phase),
            getattr(self, f"H_{band}"),
            getattr(self, f"G1_{band}"),
            getattr(self, f"G2_{band}"),
        )

        # Spin part
        geo = spin_angle(obs_.ra, obs_.dec, self.alpha0, self.delta0)
        func2 = 1 - (1 - self.R) * np.abs(geo)
        func2 = 2.5 * np.log10(func2)
        return func1 + func2

    def fit_mm(self, phase, mag, mag_err, bands, ra, dec, constrain_g1g2=True):
        """Fit a phase curve using the H, G1, G2 model from Muinonen+ 2010.

        Parameters
        ----------
        phase : list or np.ndarray
            The phase angles of observation in degree.
        mag : list or np.ndarray
            The reduced magnitudes.

        constrain_g1g2 : bool
            Constrain G1 + G2 to be less than 1.
        """
        from lmfit import minimize, Parameters

        params = lmfit.Parameters()

        for band in bands:
            params.add(f"H_{band}", value=15, min=0, max=30)
            params.add(f"G1_{band}", value=0.15, min=0, max=1.0)
            # constrain_g1g2 = False
            # params.add(f"G1_{band}", value=0.15)

            # Add delta to implement inequality constraint
            # https://lmfit.github.io/lmfit-py/constraints.html#using-inequality-constraints
            if constrain_g1g2:
                params.add(f"delta_{band}", min=0, value=0.5, max=1, vary=True)
                params.add(f"G2_{band}", expr=f"delta_{band}-G1_{band}", min=0.0, max=1)

            else:
                params.add(f"G2_{band}", value=0.15, min=0.0, max=1)
                # params.add(f"G2_{band}", value=0.15)

        params.add(f"R", value=0.8, min=0.1, max=1)
        params.add(f"alpha0", value=np.pi, min=0.0, max=2 * np.pi)
        params.add(f"delta0", value=0, min=-np.pi / 2, max=np.pi / 2)

        # params = Parameters()
        # params.add('amp', value=10)
        # params.add('decay', value=0.007)
        # params.add('phase', value=0.2)
        # params.add('frequency', value=3.0)
        #
        # Ensure that errors are well behaved
        # mag_err[mag_err == 0] = np.nanmean(mag_err)
        # mag_err = np.ones(mag.shape)
        weights = weights_from_phase(phase)  # used for weighting
        # weights = np.ones(mag.shape)

        out = minimize(
            residual,
            params,
            args=(np.radians(phase), mag, weights, bands, ra, dec),
            method="least_squares",
            jac="2-point",
            # fit_kws={
            loss="soft_l1",
            # },
        )

        for name, param in out.params.items():
            if name in ["alpha0", "delta0"]:
                setattr(self, name, np.degrees(param.value))
                # setattr(self, f"{name}_err", np.degrees(param.stderr))
            setattr(self, name, param.value)
            setattr(self, f"{name}_err", param.stderr)


def build_eqs_for_spins(x, filters=[], ph=[], ra=[], dec=[], rhs=[]):
    """Build the system of equations to solve using the HG1G2 + spin model

    Parameters
    ----------
    x: list
        List of parameters to fit for
    filters: np.array
        Array of size N containing the filtername for each measurement
    ph: np.array
        Array of size N containing phase angles
    ra: np.array
        Array of size N containing the RA (radian)
    dec: np.array
        Array of size N containing the Dec (radian)
    rhs: np.array
        Array of size N containing the actual measurements (magnitude)

    Returns
    ----------
    out: np.array
        Array of size N containing (model - y)

    Notes
    ----------
    the input `x` should start with filter independent variables,
    that is (R, alpha, delta), followed by filter dependent variables,
    that is (H, G1, G2). For example with two bands g & r:

    ```
    x = [
        R, alpha, delta,
        h_g, g_1_g, g_2_g,
        h_r, g_1_r, g_2_r
    ]
    ```

    """
    R, alpha, delta = x[0:3]
    filternames = np.unique(filters)

    params = x[3:]
    nparams = len(params) / len(filternames)
    assert int(nparams) == nparams, "You need to input all parameters for all bands"

    params_per_band = np.reshape(params, (len(filternames), int(nparams)))
    eqs = []
    for index, filtername in enumerate(filternames):
        mask = filters == filtername

        myfunc = (
            func_hg1g2_with_spin(
                np.vstack([ph[mask].tolist(), ra[mask].tolist(), dec[mask].tolist()]),
                params_per_band[index][0],
                params_per_band[index][1],
                params_per_band[index][2],
                R,
                alpha,
                delta,
            )
            - rhs[mask]
        )

        eqs = np.concatenate((eqs, myfunc))

    return np.ravel(eqs)


def func_hg1g2(ph, h, g1, g2):
    """Return f(H, G1, G2) part of the lightcurve in mag space

    Parameters
    ----------
    ph: array-like
        Phase angle in radians
    h: float
        Absolute magnitude in mag
    G1: float
        G1 parameter (no unit)
    G2: float
        G2 parameter (no unit)

    Returns
    ----------
    out: array of floats
        H - 2.5 log(f(G1G2))
    """
    from sbpy.photometry import HG1G2

    # Standard G1G2 part
    func1 = (
        g1 * HG1G2._phi1(ph) + g2 * HG1G2._phi2(ph) + (1 - g1 - g2) * HG1G2._phi3(ph)
    )
    func1 = -2.5 * np.log10(func1)

    return h + func1


def spin_angle(ra, dec, alpha0, delta0):
    return np.sin(dec) * np.sin(delta0) + np.cos(dec) * np.cos(delta0) * np.cos(
        ra - alpha0
    )


def func_hg1g2_with_spin(pha, h, g1, g2, R, alpha0, delta0):
    """Return f(H, G1, G2, R, alpha0, delta0) part of the lightcurve in mag space

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
    alpha0: float
        RA of the spin (radian)
    delta0: float
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
    func1 = func_hg1g2(ph, h, g1, g2)

    # Spin part
    geo = spin_angle(ra, dec, alpha0, delta0)
    func2 = 1 - (1 - R) * np.abs(geo)
    func2 = 2.5 * np.log10(func2)

    return func1 + func2


def residual(pars, phase, mag, weights, bands, ra, dec):
    """Build the system of equations to solve using the HG1G2 + spin model

    ```
    x = [
        R, alpha, delta,
        h_g, g_1_g, g_2_g,
        h_r, g_1_r, g_2_r
    ]
    ```

    """
    R, alpha, delta = pars["R"], pars["alpha0"], pars["delta0"]
    # filternames = np.unique(bands)

    # params = x[3:]
    params = [value for name, value in pars.items() if not name.startswith("delta")]
    params += [delta]
    params = params[:-3]
    filternames = [
        param.name.split("_")[1]
        for param in pars.values()
        if param.name.startswith("H")
    ]
    # params = list(pars.values())[:-3]
    nparams = len(params) / len(filternames)
    # assert int(nparams) == nparams, "You need to input all parameters for all bands"

    params_per_band = np.reshape(params, (len(filternames), int(nparams)))
    eqs = []
    for index, filtername in enumerate(filternames):
        mask = bands == filtername

        myfunc = (
            func_hg1g2_with_spin(
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


def weights_from_phase(phase):
    """Compute weights based on phase angle in degree"""

    # 0-2 - weight 5
    # 2-5 - weight 4
    # 5-10 - weight 3
    # 10-20 - weight 2
    # 20- weight 1
    weights = np.ones_like(phase)
    weights[phase < 2] = 5
    weights[(phase >= 2) & (phase < 5)] = 4
    weights[(phase >= 5) & (phase < 10)] = 3
    weights[(phase >= 10) & (phase < 20)] = 2
    weights[phase >= 20] = 1
    return weights
