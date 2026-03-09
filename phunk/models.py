"""phunk: Implementation of photometric models"""

import lmfit
import numpy as np
from sbpy import photometry as phot

from phunk.geometry import cos_aspect_angle, rotation_phase, subobserver_longitude
from phunk.equations import residual_shg1g2, residual_socca, func_shg1g2, func_socca
from phunk.reparametrization import (
    build_bounds,
    dict_to_lmfit,
    lmfit_to_dict,
    parameter_remapping,
    propagate_errors,
)

MODELS = ["LinExp", "HG", "HG12", "HG12S", "HG1G2", "sHG1G2", "SOCCA"]


class HG:
    """HG model from Bowell+ 1989"""

    def __init__(self, H=np.nan, G=np.nan, bands=None):
        self.NAME = "HG"
        self.PARAMS = ("H", "G")

        if bands is None:
            bands = [""]

        for band in bands:
            setattr(self, f"H{band}", H)
            setattr(self, f"G{band}", G)

    def eval(self, phase, H=None, G=None, band=""):
        """Evaluate HG model for given phase and parameters."""
        return phot.HG.evaluate(
            pha=np.radians(phase),
            hh=H if H is not None else getattr(self, f"H{band}"),
            gg=G if G is not None else getattr(self, f"G{band}"),
        )

    def is_fittable(self, pc):
        """Check if phase curve can be fit with this model."""
        _has_phase_and_mag(pc, self.NAME)
        return True

    def fit(self, pc, weights=None):
        """Fit a phase curve using the HG model."""

        def eval_fit(phase, H, G):
            """Evaluation function for fitting. Required for lmfit."""
            return phot.HG.evaluate(np.radians(phase), H, G)

        model = lmfit.Model(eval_fit)
        params = lmfit.Parameters()

        for band in pc.bands:
            params.add("H", value=15, min=-3, max=30)
            params.add("G", value=0.15, min=0, max=1.0)

            phase = pc.phase[pc.band == band]
            mag = pc.mag[pc.band == band]

            if weights is not None:
                weights_band = weights[pc.band == band]
            else:
                weights_band = np.ones(phase.shape)

            # And fit
            result = model.fit(
                mag,
                params,
                phase=phase,
                weights=weights_band,
                method="least_squares",
                fit_kws={"loss": "soft_l1"},
            )

            for param in self.PARAMS:
                setattr(self, "".join([param, band]), result.params[param].value)
                setattr(
                    self, f"{''.join([param, band])}_err", result.params[param].stderr
                )

        pc.fitted_models.add("HG")


class HG1G2:
    """HG1G2 model from Muinonen+ 2010."""

    def __init__(self, H=np.nan, G1=np.nan, G2=np.nan, bands=None):
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

        for band in pc.bands:
            params = lmfit.Parameters()
            params.add(f"H", value=15, min=-3, max=30)
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
                weights_band = weights[pc.band == band]
            else:
                weights_band = np.ones(phase.shape)

            result = model.fit(
                mag,
                params,
                phase=phase,
                weights=weights_band,
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

    def __init__(self, H=np.nan, G12=np.nan, bands=None):
        self.NAME = "HG12"
        self.PARAMS = ("H", "G12")

        if bands is None:
            bands = [""]

        for band in bands:
            setattr(self, f"H{band}", H)
            setattr(self, f"G12{band}", G12)

    def eval(self, phase, H=None, G12=None, band=""):
        """H,G1,G2 phase curve model."""
        return phot.HG12.evaluate(
            np.radians(phase),
            H if H is not None else getattr(self, f"H{band}"),
            G12 if G12 is not None else getattr(self, f"G12{band}"),
        )

    def is_fittable(self, pc):
        """Check if phase curve can be fit."""
        _has_phase_and_mag(pc, self.NAME)
        return True

    def fit(self, pc, weights=None):
        """Fit a phase curve using the H, G1, G2 model from Muinonen+ 2010.

        Parameters
        ----------
        phase : list or np.ndarray
            The phase angles of observation in degree.
        mag : list or np.ndarray
            The reduced magnitudes.
        """

        def eval_fit(phase, H, G12):
            """Evaluation function for fitting. Required for lmfit."""
            return phot.HG12.evaluate(np.radians(phase), H, G12)

        model = lmfit.Model(eval_fit)

        for band in pc.bands:
            params = lmfit.Parameters()
            params.add("H", value=15, min=-3, max=30)
            params.add("G12", value=0.3, min=0, max=1.0)

            phase = pc.phase[pc.band == band]
            mag = pc.mag[pc.band == band]

            if weights is not None:
                weights_band = weights[pc.band == band]
            else:
                weights_band = np.ones(phase.shape)

            # And fit
            result = model.fit(
                mag,
                params,
                phase=phase,
                weights=weights_band,
                method="least_squares",
                fit_kws={"loss": "soft_l1"},
            )

            for param in self.PARAMS:
                setattr(self, "".join([param, band]), result.params[param].value)
                setattr(
                    self, f"{''.join([param, band])}_err", result.params[param].stderr
                )

        pc.fitted_models.add("HG12")


class HG12S:
    """HG12* model from Penttilä+ 2016."""

    def __init__(self, H=np.nan, G12S=np.nan, bands=None):
        self.NAME = "HG12S"
        self.PARAMS = ("H", "G12S")

        if bands is None:
            bands = [""]

        for band in bands:
            setattr(self, f"H{band}", H)
            setattr(self, f"G12S{band}", G12S)

    def eval(self, phase, H=None, G12S=None, band=""):
        """H,G1,G2 phase curve model."""
        return phot.HG12_Pen16.evaluate(
            np.radians(phase),
            H if H is not None else getattr(self, f"H{band}"),
            G12S if G12S is not None else getattr(self, f"G12S{band}"),
        )

    def is_fittable(self, pc):
        """Check if phase curve can be fit."""
        _has_phase_and_mag(pc, self.NAME)
        return True

    def fit(self, pc, weights=None):
        """Fit a phase curve using the H, G1, G2 model from Muinonen+ 2010.

        Parameters
        ----------
        phase : list or np.ndarray
            The phase angles of observation in degree.
        mag : list or np.ndarray
            The reduced magnitudes.
        """

        def eval_fit(phase, H, G12S):
            """Evaluation function for fitting. Required for lmfit."""
            return phot.HG12_Pen16.evaluate(np.radians(phase), H, G12S)

        model = lmfit.Model(eval_fit)

        for band in pc.bands:
            params = lmfit.Parameters()
            params.add("H", value=15, min=0, max=30)
            params.add("G12S", value=0.3, min=0, max=1.0)

            phase = pc.phase[pc.band == band]
            mag = pc.mag[pc.band == band]

            if weights is not None:
                weights_band = weights[pc.band == band]
            else:
                weights_band = np.ones(phase.shape)

            result = model.fit(
                mag,
                params,
                phase=phase,
                weights=weights_band,
                method="least_squares",
                fit_kws={"loss": "soft_l1"},
            )

            for param in self.PARAMS:
                setattr(self, "".join([param, band]), result.params[param].value)
                setattr(
                    self, f"{''.join([param, band])}_err", result.params[param].stderr
                )

        pc.fitted_models.add("HG12S")


class LinExp:
    """Linear-Exponential Modeling from Muinonen+ 2002"""

    def __init__(self, a=np.nan, b=np.nan, d=np.nan, k=np.nan, bands=None):
        self.NAME = "LinExp"
        self.PARAMS = ("a", "b", "d", "k")

        if bands is None:
            bands = [""]

        for band in bands:
            setattr(self, f"a{band}", a)
            setattr(self, f"b{band}", b)
            setattr(self, f"d{band}", d)
            setattr(self, f"k{band}", k)

    def eval(self, phase, a=None, b=None, d=None, k=None, band=""):
        """Evaluate photometric phase curve model."""
        a = (a if a is not None else getattr(self, f"a{band}"),)
        b = (b if b is not None else getattr(self, f"b{band}"),)
        d = (d if d is not None else getattr(self, f"d{band}"),)
        k = (k if k is not None else getattr(self, f"k{band}"),)
        return a * np.exp(-phase / d) + b + k * phase

    def is_fittable(self, pc):
        """Check if phase curve can be fit."""
        _has_phase_and_mag(pc, self.NAME)
        return True

    def fit(self, pc, weights=None):
        """Fit a phase curve using the H, G1, G2 model from Muinonen+ 2010.

        Parameters
        ----------
        phase : list or np.ndarray
            The phase angles of observation in degree.
        mag : list or np.ndarray
            The reduced magnitudes.
        """

        def eval_fit(phase, a, b, d, k):
            """Evaluation function for fitting. Required for lmfit."""
            return a * np.exp(-phase / d) + b + k * phase

        model = lmfit.Model(eval_fit)

        for band in pc.bands:
            params = lmfit.Parameters()
            params.add("a", value=-15, max=30)
            params.add("b", value=-10)
            params.add("k", value=0.3)
            params.add("d", value=-3)

            phase = pc.phase[pc.band == band]
            mag = pc.mag[pc.band == band]

            if weights is not None:
                weights_band = weights[pc.band == band]
            else:
                weights_band = np.ones(phase.shape)

            result = model.fit(
                mag,
                params,
                phase=phase,
                weights=weights_band,
                method="least_squares",
                fit_kws={"loss": "soft_l1"},
            )

            for param in self.PARAMS:
                setattr(self, "".join([param, band]), result.params[param].value)
                setattr(
                    self, f"{''.join([param, band])}_err", result.params[param].stderr
                )

        pc.fitted_models.add("LinExp")


class sHG1G2:
    """sHG1G2 phase curve model."""

    def __init__(
        self,
        H=np.nan,
        G1=np.nan,
        G2=np.nan,
        alpha=np.nan,
        delta=np.nan,
        R=np.nan,
        bands=None,
    ):
        """sHG1G2 phase curve model.

        Parameters
        ----------
        H: float
            Absolute magnitude in mag
        G1: float
            G1 parameter (no unit)
        G2: float
            G2 parameter (no unit)
        alpha: float
            RA of the spin (radian)
        delta: float
            Dec of the spin (radian)
        R: float
            Oblateness (no units)
        """
        self.NAME = "sHG1G2"
        self.PARAMS = ("H", "G1", "G2", "alpha", "delta", "R")

        if bands is None:
            bands = [""]

        for band in bands:
            setattr(self, f"H{band}", H)
            setattr(self, f"G1{band}", G1)
            setattr(self, f"G2{band}", G2)

        self.alpha = alpha
        self.delta = delta
        self.R = R

    def eval(
        self,
        phase,
        ra,
        dec,
        band="",
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

        H = self.H if not band else getattr(self, f"H{band}")
        G1 = self.G1 if not band else getattr(self, f"G1{band}")
        G2 = self.G2 if not band else getattr(self, f"G2{band}")

        alpha = self.alpha if alpha is None else alpha
        delta = self.delta if delta is None else delta
        R = self.R if R is None else R

        func = func_shg1g2([phase, ra, dec], H, G1, G2, R, alpha, delta)

        return func

    def is_fittable(self, pc):
        """Check whether phase curve can be fit with sHG1G2 model."""
        _has_phase_and_mag(pc, self.NAME)

        # Do we have RA and Dec?
        if pc.ra is None or pc.dec is None:
            assert pc.epoch is not None and pc.target is not None, (
                f"Neither RA/Dec nor target nor observation epoch provided. Cannot fit {self.NAME} model."
            )

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

        # Fitting all bands at once with sHG1G2
        params.add(f"R", value=0.8, min=0.1, max=1)
        params.add(f"alpha", value=np.pi, min=0.0, max=2 * np.pi)
        params.add(f"delta", value=0, min=-np.pi / 2, max=np.pi / 2)

        if weights is None:
            weights = np.ones(pc.mag.shape)

        out = lmfit.minimize(
            residual_shg1g2,
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
            loss="soft_l1",
        )

        for name, param in out.params.items():
            if name in ["alpha", "delta"]:
                setattr(self, name, np.degrees(param.value))
                continue
            setattr(self, name, param.value)
            setattr(self, f"{name}_err", param.stderr)
        pc.fitted_models.add("sHG1G2")


class SOCCA:
    """SOCCA phase curve model."""

    def __init__(
        self,
        H=np.nan,
        G1=np.nan,
        G2=np.nan,
        period=np.nan,
        alpha=np.nan,
        delta=np.nan,
        a_b=np.nan,
        a_c=np.nan,
        W0=np.nan,
        t0=np.nan,
        bands=None,
        p0=None,
        remap=False,
    ):
        """SOCCA phase curve model.

        Parameters
        ----------
        H: float
            Absolute magnitude in mag
        G1: float
            G1 parameter (no unit)
        G2: float
            G2 parameter (no unit)
        period: float
            Sidereal rotation period (days)
        alpha: float
            RA of the spin (degree)
        delta: float
            Dec of the spin (degree)
        a_b: float
            Triaxial axes a/b ratio (no units)
        a_c: float
            Triaxial axes a/c ratio (no units)
        W0: float
            Initial rotation angle (degree)
        t0: float
            Reference epoch for rotation (JD)
        p0: dict
            Dictionary with guess estimates
        """
        self.NAME = "SOCCA"
        self.PARAMS = (
            "H",
            "G1",
            "G2",
            "period",
            "alpha",
            "delta",
            "a_b",
            "a_c",
            "W0",
            "t0",
        )

        if bands is None:
            bands = [""]

        for band in bands:
            setattr(self, f"H{band}", H)
            setattr(self, f"G1{band}", G1)
            setattr(self, f"G2{band}", G2)

        self.period = period
        self.alpha = alpha
        self.delta = delta
        self.a_b = a_b
        self.a_c = a_c
        self.W0 = W0
        self.t0 = t0

        self.p0 = p0
        if isinstance(p0, dict):
            for k, v in p0.items():
                setattr(self, k, v)
        self.remap = remap

    def eval(
        self,
        phase,
        ra,
        dec,
        epoch,
        ra_s,
        dec_s,
        band="",
        H=None,
        G1=None,
        G2=None,
        period=None,
        alpha=None,
        delta=None,
        a_b=None,
        a_c=None,
        W0=None,
        t0=None,
    ):
        """SOCCA phase curve model.

        Return f(H, G1, G2, period, alpha, delta, a_b, a_c, W0, t0) part of the lightcurve in mag space

        Parameters
        ----------
        H: float
            Absolute magnitude in mag
        G1: float
            G1 parameter (no unit)
        G2: float
            G2 parameter (no unit)
        period: float
            Sidereal rotation period (days)
        alpha: float
            Right ascension of the spin axis (degree)
        delta: float
            Declination of the spin axis (degree)
        a_b: float
            Triaxial axes a/b ratio (no units)
        a_c: float
            Triaxial axes a/c ratio (no units)
        W0: float
            Initial rotation angle (degree)
        t0: float
            Reference epoch for rotation (JD)

        Returns
        ----------
        out: array of floats
            H - 2.5 log(f(G1G2)) - 2.5 log(f(R, spin))
        """
        phase = np.radians(phase)
        ra = np.radians(ra)
        dec = np.radians(dec)

        ra_s = np.radians(ra_s)
        dec_s = np.radians(dec_s)

        H = self.H if not band else getattr(self, f"H{band}")
        G1 = self.G1 if not band else getattr(self, f"G1{band}")
        G2 = self.G2 if not band else getattr(self, f"G2{band}")

        period = self.period if period is None else period
        alpha = np.radians(self.alpha) if alpha is None else np.radians(alpha)
        delta = np.radians(self.delta) if delta is None else np.radians(delta)
        a_b = self.a_b if a_b is None else a_b
        a_c = self.a_c if a_c is None else a_c
        W0 = np.radians(self.W0) if W0 is None else np.radians(W0)
        t0 = self.t0 if t0 is None else t0

        func = func_socca(
            [phase, ra, dec, epoch, ra_s, dec_s],
            H,
            G1,
            G2,
            alpha,
            delta,
            period,
            a_b,
            a_c,
            W0,
        )

        return func

    def is_fittable(self, pc):
        """Check whether phase curve can be fit with SOCCA model."""
        _has_phase_and_mag(pc, self.NAME)

        # Do we have RA and Dec?
        if pc.ra is None or pc.dec is None:
            assert pc.epoch is not None and pc.target is not None, (
                f"Neither RA/Dec nor target nor observation epoch provided. Cannot fit {self.NAME} model."
            )

            print("No RA or Dec provided but epoch and target found.")
            pc.get_ephems()

        return True

    def fit(self, pc, weights=None):
        """Fit a phase curve using the SOCCA model."""

        ### Initialize with unbounded parameters ###
        params = lmfit.Parameters()
        for band in set(pc.band):
            params.add(f"H{band}", value=self.H)
            params.add(f"G1{band}", value=self.G1)
            params.add(f"G2{band}", value=self.G2)
            
        # Fitting all bands at once with SOCCA
        params.add(
            "alpha",
            value=np.radians(self.alpha),
        )
        params.add(
            "delta",
            value=np.radians(self.delta),
        )
        params.add("period", value=self.period)
        params.add("a_b", value=1.15)
        params.add("a_c", value=1.5)
        params.add("W0", value=np.radians(self.W0))
        params.add("t0", value=self.t0, vary=False)

        ### Reparametrize part ###
        if self.remap:
            params_dict = lmfit_to_dict(params)
            latent_dict = parameter_remapping(params_dict, physical_to_latent=True)
            params = dict_to_lmfit(latent_dict, params)

        lower, upper = build_bounds(remap=self.remap)
        if self.remap:
            fbase = ["H", "u_G1", "u_G2"]
        else:
            fbase = ["H", "G1", "G2"]
        for band in pc.bands:
            for base in fbase:
                key = f"{base}{band}"
                lower[key] = lower[base]
                upper[key] = upper[base]
        for base in fbase:
            lower.pop(base, None)
            upper.pop(base, None)

        for name in params.keys():
            params[name].min = lower[name]
            params[name].max = upper[name]
        if weights is None:
            weights = np.ones(pc.mag.shape)

        out = lmfit.minimize(
            residual_socca,
            params,
            args=(
                np.radians(pc.phase),
                pc.mag,
                weights,
                pc.band,
                np.radians(pc.ra),
                np.radians(pc.dec),
                pc.epoch,
                np.radians(pc.ra_s),
                np.radians(pc.dec_s),
                self.remap,
            ),
            method="least_squares",
            jac="2-point",
            loss="soft_l1",
        )

        out_min = out.params

        if self.remap:
            latent_min = lmfit_to_dict(out_min)
            physical_min = parameter_remapping(latent_min, physical_to_latent=False)
        else:
            physical_min = lmfit_to_dict(out_min)

        err_dict = {lname: lparam.stderr for lname, lparam in out.params.items()}
        if self.remap:
            err_dict = propagate_errors(latent_min, err_dict, filt_names=set(pc.bands))
        for name, param in physical_min.items():
            if name in ["alpha", "delta", "W0"]:
                setattr(self, name, np.degrees(param))
                setattr(self, f"{name}_err", np.degrees(err_dict[name]))
                continue
            setattr(self, name, param)
            setattr(self, f"{name}_err", err_dict[name])

        pc.fitted_models.add("SOCCA")


def _has_phase_and_mag(pc, model):
    """Check that phase curve as magnitude and phase attributes."""
    # We need at least magnitudes
    assert pc.mag is not None, f"No magnitudes provided. Cannot fit {model} model."

    if pc.phase is None:
        assert pc.epoch is not None and pc.target is not None, (
            f"Neither phase nor target nor observation epoch provided. Cannot fit {model} model."
        )

        print("No phase angles provided but epoch and target found.")
        pc.get_ephems()

    # If we have matching phase angles, all is good
    if pc.phase is not None:
        assert len(pc.mag) == len(pc.phase), (
            f"Provided 'phase' and 'mag' have different lengths: {len(pc.phase)} vs {len(pc.mag)}. "
            f"Cannot fit {model} model."
        )
        return True
