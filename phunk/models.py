"""phunk: Implementation of photometric models"""

import lmfit
import numpy as np
from sbpy import photometry as phot

MODELS = ["LinExp", "HG", "HG12", "HG12S", "HG1G2", "sHG1G2"]


class HG:
    """HG model from Bowell+ 1989"""

    def __init__(self, H=10, G=0.15, bands=None):
        self.NAME = "HG"
        self.PARAMS = ("H", "G")

        self.H = H
        self.G = G

    def eval(self, phase, H=None, G=None):
        """Evaluate HG model for given phase and parameters."""
        return phot.HG.evaluate(
            pha=np.radians(phase),
            hh=H if H is not None else self.H,
            gg=G if G is not None else self.G,
        )

    def is_fittable(self, pc):
        """Check if phase curve can be fit with this model."""
        _has_phase_and_mag(pc, self.NAME)
        return True

    def fit(self, pc):
        """Fit a phase curve using the HG model."""

        # The model function
        model = lmfit.Model(self.eval)

        # The fit parameters
        params = lmfit.Parameters()
        params.add("H", value=15, min=0, max=30)
        params.add("G", value=0.15, min=0, max=1.0)

        # And fit
        result = model.fit(
            pc.mag,
            params,
            phase=pc.phase,
            method="least_squares",
            fit_kws={
                "loss": "soft_l1",
            },
            # weights=weights,
        )

        for param in self.PARAMS:
            setattr(self, param, result.params[param].value)
            setattr(self, f"{param}_err", result.params[param].stderr)

        pc.fitted_models.add("HG")


class HG1G2:
    """HG1G2 model from Muinonen+ 2010."""

    def __init__(self, H=15.0, G1=0.15, G2=0.15, bands=None):
        self.NAME = "HG1G2"
        self.PARAMS = ("H", "G1", "G2")

        if bands is None:
            bands = [""]

        for band in bands:
            setattr(self, f"H{band}", H)
            setattr(self, f"G1{band}", G1)
            setattr(self, f"G2{band}", G2)

    def eval(self, phase, H=None, G1=None, G2=None, band=""):
        """H,G1,G2 phase curve model."""
        return phot.HG1G2.evaluate(
            np.radians(phase),
            H if H is not None else getattr(self, f"H{band}"),
            G1 if G1 is not None else getattr(self, f"G1{band}"),
            G2 if G2 is not None else getattr(self, f"G2{band}"),
        )

    def is_fittable(self, pc):
        """Check if phase curve can be fit."""
        _has_phase_and_mag(pc, self.NAME)
        return True

    def fit(self, pc, weights=None, constrain_g1g2=False):
        """Fit a phase curve using the H, G1, G2 model from Muinonen+ 2010.

        Parameters
        ----------
        """

        def eval_fit(phase, H, G1, G2):
            """Evaluation function for fitting. Required for lmfit."""
            return phot.HG1G2.evaluate(np.radians(phase), H, G1, G2)

        model = lmfit.Model(eval_fit)

        for band in set(pc.band):
            params = lmfit.Parameters()
            params.add(f"H", value=15, min=0, max=30)
            params.add(f"G1", value=0.15, min=0, max=1.0)

            if constrain_g1g2:
                # Add delta to implement inequality constraint
                # https://lmfit.github.io/lmfit-py/constraints.html#using-inequality-constraints
                params.add("delta", min=0, value=0.5, max=1, vary=True)
                params.add(f"G2", expr="delta-G1", min=0.0, max=1)
            else:
                params.add(f"G2", value=0.15, min=0, max=1.0)

            phase = pc.phase[pc.band == band]
            mag = pc.mag[pc.band == band]

            if weights is not None:
                weights = weights[pc.band == band]

            result = model.fit(
                mag,
                params,
                phase=phase,
                weights=weights,
                method="least_squares",
                fit_kws={"loss": "soft_l1"},
            )

            for param in self.PARAMS:
                setattr(self, "".join([param, band]), result.params[param].value)
                setattr(
                    self, f"{''.join([param, band])}_err", result.params[param].stderr
                )

        pc.fitted_models.add("HG1G2")


class HG12:
    """HG12 model from Muinonen+ 2010."""

    def __init__(self, H=15.0, G12=0.15, bands=None):
        self.NAME = "HG12"
        self.PARAMS = ("H", "G12")

        self.H = H
        self.G1 = G12

    def eval(self, phase, H=None, G12=None):
        """H,G1,G2 phase curve model."""
        phase = np.radians(phase)

        if H is None:
            H = self.H
        if G12 is None:
            G12 = self.G12

        return phot.HG12.evaluate(phase, H, G12)

    def is_fittable(self, pc):
        """Check if phase curve can be fit."""
        _has_phase_and_mag(pc, self.NAME)
        return True

    def fit(self, pc):
        """Fit a phase curve using the H, G1, G2 model from Muinonen+ 2010.

        Parameters
        ----------
        phase : list or np.ndarray
            The phase angles of observation in degree.
        mag : list or np.ndarray
            The reduced magnitudes.
        """

        # The model function
        model = lmfit.Model(self.eval)

        # The fit parameters
        params = lmfit.Parameters()
        params.add("H", value=15, min=0, max=30)
        params.add("G12", value=0.3, min=0, max=1.0)

        # And fit
        result = model.fit(
            pc.mag,
            params,
            phase=pc.phase,
            method="least_squares",
            fit_kws={
                "loss": "soft_l1",
            },
        )

        for param in self.PARAMS:
            setattr(self, param, result.params[param].value)
            setattr(self, f"{param}_err", result.params[param].stderr)

        pc.fitted_models.add("HG12")


class HG12S:
    """HG12* model from Penttilä+ 2016."""

    def __init__(self, H=15.0, G12S=0.15, bands=None):
        self.NAME = "HG12S"
        self.PARAMS = ("H", "G12S")

        self.H = H
        self.G12S = G12S

    def eval(self, phase, H=None, G12S=None):
        """H,G1,G2 phase curve model."""
        phase = np.radians(phase)

        if H is None:
            H = self.H
        if G12S is None:
            G12S = self.G12S

        return phot.HG12_Pen16.evaluate(phase, H, G12S)

    def is_fittable(self, pc):
        """Check if phase curve can be fit."""
        _has_phase_and_mag(pc, self.NAME)
        return True

    def fit(self, pc):
        """Fit a phase curve using the H, G1, G2 model from Muinonen+ 2010.

        Parameters
        ----------
        phase : list or np.ndarray
            The phase angles of observation in degree.
        mag : list or np.ndarray
            The reduced magnitudes.
        """

        # The model function
        model = lmfit.Model(self.eval)

        # The fit parameters
        params = lmfit.Parameters()
        params.add("H", value=15, min=0, max=30)
        params.add("G12S", value=0.3, min=0, max=1.0)

        # And fit
        result = model.fit(
            pc.mag,
            params,
            phase=pc.phase,
            method="least_squares",
            fit_kws={
                "loss": "soft_l1",
            },
        )

        for param in self.PARAMS:
            setattr(self, param, result.params[param].value)
            setattr(self, f"{param}_err", result.params[param].stderr)

        pc.fitted_models.add("HG12S")


class LinExp:
    """Linear-Exponential Modeling from Muinonen+ 2002"""

    def __init__(self, a=15.0, b=0.15, d=0.15, k=1, bands=None):
        self.NAME = "LinExp"
        self.PARAMS = ("a", "b", "d", "k")

        self.a = a
        self.b = b
        self.d = d
        self.k = k

    def eval(self, phase, a=None, b=None, d=None, k=None):
        """Evaluate photometric phase curve model."""

        if a is None:
            a = self.a
        if b is None:
            b = self.b
        if d is None:
            d = self.d
        if k is None:
            k = self.k

        return a * np.exp(-phase / d) + b + k * phase

    def is_fittable(self, pc):
        """Check if phase curve can be fit."""
        _has_phase_and_mag(pc, self.NAME)
        return True

    def fit(self, pc):
        """Fit a phase curve using the H, G1, G2 model from Muinonen+ 2010.

        Parameters
        ----------
        phase : list or np.ndarray
            The phase angles of observation in degree.
        mag : list or np.ndarray
            The reduced magnitudes.
        """

        # The model function
        model = lmfit.Model(self.eval)

        # The fit parameters
        params = lmfit.Parameters()
        params.add("a", value=-15, max=30)
        params.add("b", value=-10)
        params.add("k", value=0.3)
        params.add("d", value=-3)

        # And fit
        result = model.fit(
            pc.mag,
            params,
            phase=pc.phase,
            method="least_squares",
            fit_kws={
                "loss": "soft_l1",
            },
        )

        for param in self.PARAMS:
            setattr(self, param, result.params[param].value)
            setattr(self, f"{param}_err", result.params[param].stderr)

        pc.fitted_models.add("LinExp")


class sHG1G2:
    """sHG1G2 phase curve model."""

    def __init__(self, bands=None):
        """
        Provide ra and dec in degrees
        """
        self.NAME = "sHG1G2"
        self.PARAMS = ("H", "G1", "G2", "alpha", "delta", "R")

        # self.alpha = alpha
        # self.delta = delta
        # self.R = R

    def eval(
        self,
        phase,
        ra,
        dec,
        band,
        H=None,
        G1=None,
        G2=None,
        alpha=None,
        delta=None,
        R=None,
    ):
        """sHG1G2 phase curve model. Adapted implementation from J. Peloton.

        Return f(H, G1, G2, R, alpha, delta) part of the lightcurve in mag space

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
        alpha: float
            RA of the spin (radian)
        delta: float
            Dec of the spin (radian)

        Returns
        ----------
        out: array of floats
            H - 2.5 log(f(G1G2)) - 2.5 log(f(R, spin))
        """
        phase = np.radians(phase)
        ra = np.radians(ra)
        dec = np.radians(dec)

        if H is None:
            H = self.H
        if G1 is None:
            G1 = self.G1
        if G2 is None:
            G2 = self.G2
        if alpha is None:
            alpha = self.alpha
        if delta is None:
            delta = self.delta
        if R is None:
            R = self.R

        # Standard HG1G2 part: h + f(alpha, G1, G2)
        func1 = phot.HG1G2().evaluate(
            phase,
            getattr(self, f"H{band}"),
            getattr(self, f"G1{band}"),
            getattr(self, f"G2{band}"),
        )

        # Spin part
        geo = spin_angle(ra, dec, alpha, delta)
        func2 = 1 - (1 - self.R) * np.abs(geo)
        func2 = 2.5 * np.log10(func2)
        return func1 + func2

    def is_fittable(self, pc):
        """Check whether phase curve can be fit with sHG1G2 model."""
        _has_phase_and_mag(pc, self.NAME)

        # Do we have RA and Dec?
        if pc.ra is None or pc.dec is None:
            assert (
                pc.epoch is not None and pc.target is not None
            ), f"Neither RA/Dec nor target nor observation epoch provided. Cannot fit {self.NAME} model."

            print("No RA or Dec provided but epoch and target found.")
            pc.get_ephems()

        return True

    def fit(self, pc, weights=None, constrain_g1g2=False):
        """Fit a phase curve using the sHG1G2 model."""

        params = lmfit.Parameters()

        for band in set(pc.band):
            params.add(f"H{band}", value=15, min=0, max=30)
            params.add(f"G1{band}", value=0.15, min=0, max=1.0)

            # Add delta to implement inequality constraint
            # https://lmfit.github.io/lmfit-py/constraints.html#using-inequality-constraints
            if constrain_g1g2:
                params.add(f"delta{band}", min=0, value=0.5, max=1, vary=True)
                params.add(f"G2{band}", expr=f"delta{band}-G1{band}", min=0.0, max=1)

            else:
                params.add(f"G2{band}", value=0.15, min=0.0, max=1)

        params.add(f"R", value=0.8, min=0.1, max=1)
        params.add(f"alpha", value=np.pi, min=0.0, max=2 * np.pi)
        params.add(f"delta", value=0, min=-np.pi / 2, max=np.pi / 2)

        if weights is None:
            weights = np.ones(pc.mag.shape)

        out = lmfit.minimize(
            residual,
            params,
            args=(
                np.radians(pc.phase),
                pc.mag,
                weights,
                pc.band,
                np.radians(pc.ra),
                np.radians(pc.dec),
            ),
            method="least_squares",
            jac="2-point",
            # fit_kws={
            loss="soft_l1",
            # },
        )

        for name, param in out.params.items():
            if name in ["alpha", "delta"]:
                setattr(self, name, np.degrees(param.value))
                continue
            setattr(self, name, param.value)
            setattr(self, f"{name}_err", param.stderr)
        pc.fitted_models.add("sHG1G2")


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
    # Standard G1G2 part
    func1 = (
        g1 * phot.HG1G2._phi1(ph)
        + g2 * phot.HG1G2._phi2(ph)
        + (1 - g1 - g2) * phot.HG1G2._phi3(ph)
    )
    func1 = -2.5 * np.log10(func1)

    return h + func1


def spin_angle(ra, dec, alpha, delta):
    return np.sin(dec) * np.sin(delta) + np.cos(dec) * np.cos(delta) * np.cos(
        ra - alpha
    )


def func_hg1g2_with_spin(pha, h, g1, g2, R, alpha, delta):
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
    func1 = func_hg1g2(ph, h, g1, g2)

    # Spin part
    geo = spin_angle(ra, dec, alpha, delta)
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


def _has_phase_and_mag(pc, model):
    """Check that phase curve as magnitude and phase attributes."""
    # We need at least magnitudes
    assert pc.mag is not None, f"No magnitudes provided. Cannot fit {model} model."

    if pc.phase is None:
        assert (
            pc.epoch is not None and pc.target is not None
        ), f"Neither phase nor target nor observation epoch provided. Cannot fit {model} model."

        print("No phase angles provided but epoch and target found.")
        pc.get_ephems()

    # If we have matching phase angles, all is good
    if pc.phase is not None:
        assert len(pc.mag) == len(pc.phase), (
            f"Provided 'phase' and 'mag' have different lengths: {len(pc.phase)} vs {len(pc.mag)}. "
            f"Cannot fit {model} model."
        )
        return True
