import numpy as np
import lmfit

# # FIXME
# # Source - https://stackoverflow.com/a/30368735
# # Posted by niekas, modified by community. See post 'Timeline' for change history
# # Retrieved 2026-03-10, License - CC BY-SA 4.0

# import warnings
# warnings.filterwarnings("error")
# # FIXME


def sigmoid(x):
    """
    Compute the sigmoid function.
    Maps any real number to the interval (0, 1).
    """
    # try:
    #     return 1 / (1 + np.exp(-x))
    # except RuntimeWarning:
    #     print(x)
    return 1 / (1 + np.exp(-x))
def sc_sigmoid(x, C=-0.429, R=1.429 + np.abs(-0.429), k=1, In=0):
    """
    Compute the scaled sigmoid function.
    Maps any real number to the interval (C, R-|C|).
    """
    return C + R / (1 + np.exp(-k * (x - In)))


def logit(x):
    """
    Compute the logit (inverse sigmoid) function.
    Maps a value in (0, 1) to the real line.
    """
    return np.log(x / (1 - x))


def sc_logit(y, C=-0.429, R=1.429 + np.abs(-0.429)):
    """
    Compute the scaled logit (inverse scaled sigmoid) function.
    Maps a value in (C, R-|C|) to the real line.
    """
    p = (y - C) / R
    return np.log(p / (1 - p))


GMIN = -0.429
GMAX = 1.429

a1, b1 = -3.9038, -0.2445
a2, b2 = -0.9635, 1.0157
a3, b3 = -0.5330, 0.0027


def compute_LU_bounds(g1):
    """
    Compute allowed interval for G2 given G1

    Parameters
    ----------
    g1: np.array
        G1 phase parameter values

    Returns
    -------
    L: np.array
        Lower bounds of G2
    U: np.array
        Upper bounds of G2
    """
    lower1 = a1 * g1 + b1
    lower2 = a3 * g1 + b3

    L = np.maximum(np.maximum(GMIN, lower1), lower2)
    U = np.minimum(GMAX, a2 * g1 + b2)

    return L, U


def build_bounds(bounds=None, remap=False):
    """
    Build lower and upper bounds for parameters with optional reparametrization.

    Parameters that are reparametrized are set to (-inf, +inf), otherwise
    default physical bounds are used.

    Order of parameters:
        H, G1, G2, alpha0, delta0, period, a_b, a_c, phi0

    Parameters
    ----------
    bounds : tuple of lists, optional
        Physical bounds ((lower list, upper list)) for each parameter.
    use_angles : bool
        If True, set bounds for spin axis coords (alpha0, delta0) to (-inf, +inf).
    use_shape : bool
        If True, set bounds for shape parameters (a_b, a_c) to (-inf, +inf).
    use_phase : bool
        If True, set bounds for phi0 to (-inf, +inf).
    use_filter_dependent : bool
        If True, set bounds for filter dependent parameters (H, G1, G2) to (-inf, +inf).

    Returns
    -------
    lower_bounds : np.ndarray
        Lower bounds for all parameters.
    upper_bounds : np.ndarray
        Upper bounds for all parameters.
    """
    physical_keys = [
        "H",
        "G1",
        "G2",
        "alpha",
        "delta",
        "period",
        "a_b",
        "a_c",
        "W0",
        "t0",
    ]
    latent_keys = [
        "H",
        "u_G1",
        "u_G2",
        "period",
        "X",
        "Y",
        "Z",
        "u_a_b",
        "u_a_c",
        "u_W0",
        "t0",
    ]

    if bounds is None:
        bounds = (
            [-3, GMIN, GMIN, 0, -np.pi / 2, 2.2 / 24.0, 1, 1, -np.pi / 2, -np.inf],
            [30, GMAX, GMAX, 2 * np.pi, np.pi / 2, 1000, 5, 5, np.pi / 2, np.inf],
        )
        lower_bounds = np.array(bounds[0])
        upper_bounds = np.array(bounds[1])
    if remap:
        bounds = (
            [
                -3,
                GMIN,
                GMIN,
                2.2 / 24.0,
                -np.inf,
                -np.inf,
                -np.inf,
                1,
                1,
                -np.pi / 2,
                -np.inf,
            ],
            [30, GMAX, GMAX, 1000, np.inf, np.inf, np.inf, 5, 5, np.pi / 2, np.inf],
        )
        lower_bounds = np.array(bounds[0])
        upper_bounds = np.array(bounds[1])

        lower_bounds[1:3] = -np.inf  # H. G1. G2
        upper_bounds[1:3] = np.inf
        # P. X. Y. Z
        lower_bounds[7:9] = -np.inf  # a_b. a_c
        upper_bounds[7:9] = np.inf

        lower_bounds[9] = -np.inf  # W0
        upper_bounds[9] = np.inf
    if remap:
        lower_dict = {k: v for k, v in zip(latent_keys, lower_bounds)}
        upper_dict = {k: v for k, v in zip(latent_keys, upper_bounds)}
    else:
        lower_dict = {k: v for k, v in zip(physical_keys, lower_bounds)}
        upper_dict = {k: v for k, v in zip(physical_keys, upper_bounds)}

    return lower_dict, upper_dict


def prop_angle_error(X, Y, Z, err_X, err_Y, err_Z):
    """
    Propagate Cartesian coordinate uncertainties to angular uncertainties.

    Computes the propagated errors on the angular spin axis coordinates
    (alpha0, delta0) from uncertainties on their Cartesian components (X, Y, Z).

    Parameters
    ----------
    X, Y, Z : float
        Cartesian coordinates of the vector.
    err_X, err_Y, err_Z : float
        1-sigma uncertainties of X, Y and Z.

    Returns
    -------
    err_alpha0 : float
        Propagated 1-sigma uncertainty of alpha0.
    err_delta0 : float or ndarray
        Propagated 1-sigma uncertainty of delta0.
    """
    dfdx = -(X * Z) / (
        np.sqrt(1 - Z**2 / (X**2 + Y**2 + Z**2)) * (X**2 + Y**2 + Z**2) ** (3 / 2)
    )
    dfdy = -(Y * Z) / (
        np.sqrt(1 - Z**2 / (X**2 + Y**2 + Z**2)) * (X**2 + Y**2 + Z**2) ** (3 / 2)
    )
    dfdz = (
        1 / np.sqrt(X**2 + Y**2 + Z**2) - Z**2 / (X**2 + Y**2 + Z**2) ** (3 / 2)
    ) * (1 / np.sqrt(1 - Z**2 / (X**2 + Y**2 + Z**2)))

    term1 = (dfdx * err_X) ** 2 + (dfdy * err_Y) ** 2 + (dfdz * err_Z) ** 2
    term2 = 2 * (
        dfdx * dfdy * err_X * err_Y
        + dfdx * dfdz * err_X * err_Z
        + dfdy * dfdz * err_Y * err_Z
    )
    err_delta0 = np.sqrt(term1 + term2)

    err_alpha0 = np.sqrt(
        (X / (X**2 + Y**2) * err_Y) ** 2
        + (Y / (X**2 + Y**2) * err_X) ** 2
        - (X * Y) / (X**2 + Y**2) ** 2 * err_Y * err_X
    )
    return err_alpha0, err_delta0


def prop_phase_error(u_phi0, err_u_phi0):
    """
    Propagates the uncertainty on u_phi0 to the corresponding
    uncertainty on the initial roation phase phi0.

    Parameters
    ----------
    u_phi0 : float
        Unconstrained initial phase
    err_u_phi0 : float
        1-sigma uncertainty on u_phi0

    Returns
    -------
    err_phi0 : float
        Propagated 1-sigma uncertainty on phi0.
    """
    err_phi0 = np.pi * sigmoid(u_phi0) * (1 - sigmoid(u_phi0)) * err_u_phi0
    return err_phi0


def prop_G1_err(u_G1, err_u_G1, R=1.429 + np.abs(-0.429)):
    """
    Propagate uncertainty from u_G1 to the G1 parameter.

    Parameters
    ----------
    u_G1 : float
        Unconstrained parameter mapped to G1.
    err_u_G1 : float
        1-sigma uncertainty on u_G1.

    Returns
    -------
    err_G1 : float
        Propagated 1-sigma uncertainty on G1.
    """
    err_G1 = R * np.exp(-u_G1) / (np.exp(-u_G1) + 1) ** 2 * err_u_G1
    return err_G1


def prop_G2_err(G1, u_G2, err_u_G2):
    """
    Propagate uncertainty to the G2 parameter.

    Parameters
    ----------
    G1 : float
        G1 phase parameter.
    u_G2 : float
        Unconstrained parameter mapped to G2.
    err_u_G2 : float
        1-sigma uncertainty on u_G2.
    err_G1 : float
        1-sigma uncertainty on G1.

    Returns
    -------
    err_G2 : float
        Propagated 1-sigma uncertainty on G2.
    """
    L, U = compute_LU_bounds(G1)
    err_G2 = (U - L) * (1 - sigmoid(u_G2)) * sigmoid(u_G2) * err_u_G2
    return err_G2


def prop_ab_err(u_a_b, err_u_a_b):
    """
    Propagate uncertainty to the a/b shape parameter.

    Parameters
    ----------
    a_b : float
        a/b shape parameter.
    u_a_b : float
        Unconstrained parameter mapped to a/b.
    Returns
    -------
    err_a_b : float
        Propagated 1-sigma uncertainty on a/b.
    """
    err_a_b = 4 * sigmoid(u_a_b) * (1 - sigmoid(u_a_b)) * err_u_a_b

    return err_a_b


def prop_ac_err(a_b, u_a_c, err_u_a_c, err_a_b):
    """
    Propagate the 1-sigma uncertainty to the a/c shape parameter.

    Parameters
    ----------
    a_b : float
        Physical a/b shape parameter.
    u_a_c : float
        Unconstrained parameter mapped to a/c.
    err_u_a_c : float
        1-sigma uncertainty on u_a_c.
    err_a_b : float
        1-sigma uncertainty on a/b.

    Returns
    -------
    err_a_c : float
        Propagated 1-sigma uncertainty on a/c.
    """
    term1 = ((1 - sigmoid(u_a_c)) * err_a_b) ** 2
    term2 = ((5 - a_b) * sigmoid(u_a_c) * (1 - sigmoid(u_a_c)) * err_u_a_c) ** 2
    term3 = (
        2
        * (1 - sigmoid(u_a_c))
        * (5 - a_b)
        * sigmoid(u_a_c)
        * (1 - sigmoid(u_a_c))
        * err_u_a_c
        * err_a_b
    )
    err_a_c = np.sqrt(term1 + term2 + term3)
    return err_a_c

def propagate_errors(params_dict, err_dict, filt_names=None):
    """
    Propagate errors from a dictionary of fit parameters to physical parameters,
    using parameter keys instead of relying on array indices.

    Parameters
    ----------
    params_dict : dict
        Dictionary of fitted parameter values, keyed by parameter name (latent or physical)
    err_dict : dict
        Dictionary of parameter errors (standard deviations), keyed by same names
    filt_names : list of str, optional
        Names of the filters, e.g., ['g', 'r', 'i'] to label H, G1, G2

    Returns
    -------
    dict
        Dictionary of propagated errors keyed by physical parameter names
    """
    if filt_names is None:
        filt_names = [""]

    out_errs = {}

    # Spin axis
    X = params_dict.get("X", None)
    Y = params_dict.get("Y", None)
    Z = params_dict.get("Z", None)
    err_X = err_dict.get("X", None)
    err_Y = err_dict.get("Y", None)
    err_Z = err_dict.get("Z", None)
    err_alpha0, err_delta0 = prop_angle_error(X, Y, Z, err_X, err_Y, err_Z)
    out_errs["alpha"] = err_alpha0
    out_errs["delta"] = err_delta0
    out_errs["period"] = err_dict.get("period", None)

    # Shape
    u_a_b = params_dict.get("u_a_b", None)
    u_a_c = params_dict.get("u_a_c", None)
    err_u_a_b = err_dict.get("u_a_b", None)
    err_u_a_c = err_dict.get("u_a_c", None)

    a_b_phys = 4 * sigmoid(u_a_b) + 1
    err_ab_phys = prop_ab_err(u_a_b=u_a_b, err_u_a_b=err_u_a_b)
    err_ac_phys = prop_ac_err(a_b=a_b_phys, u_a_c=u_a_c, err_u_a_c=err_u_a_c, err_a_b=err_u_a_b)
    out_errs["a_b"] = err_ab_phys
    out_errs["a_c"] = err_ac_phys

    # Phase
    u_phi0 = params_dict.get("u_W0", None)
    err_phi0 = err_dict.get("u_W0", None)
    err_W0 = prop_phase_error(u_phi0=u_phi0, err_u_phi0=err_phi0)
    out_errs["W0"] = err_W0

    # Filter-dependent
    for band in filt_names:
        # H = params_dict.get(f"H{band}", None)
        G1 = params_dict.get(f"u_G1{band}", None) 
        G2 = params_dict.get(f"u_G2{band}", None) 

        err_H = err_dict.get(f"H{band}", None)
        err_G1 = prop_G1_err(u_G1=G1, err_u_G1=err_dict.get(f"u_G1{band}", None))
        G1_phys = sc_sigmoid(G1)
        err_G2 = prop_G2_err(G1=G1_phys, u_G2=G2, err_u_G2=err_dict.get(f"u_G2{band}", None))

        out_errs[f"H{band}"] = err_H
        out_errs[f"G1{band}"] = err_G1
        out_errs[f"G2{band}"] = err_G2

    # t0
    out_errs["t0"] = err_dict.get("t0", None)

    return out_errs

def parameter_remapping(
    pars,
    physical_to_latent=True,
    use_angles=True,
    use_shape=True,
    use_phase=True,
    use_filter_dependent=True,
):
    """
    Convert between physical and latent parameter representations, using dictionaries.

    Parameters
    ----------
    pars : dict
        Keys are parameter names (H_g, G1_g, G2_g, alpha, delta, a_b, a_c, W0, period, t0)
    physical_to_latent : bool
        Direction of conversion. True: physical -> latent, False: latent -> physical
    use_angles, use_shape, use_phase, use_filter_dependent : bool
        Which blocks to convert

    Returns
    -------
    dict
        Dictionary of remapped parameters with the same keys
    """
    pars_out = {}
    pars = dict(pars)
    
    # -----------------------------
    # Angles / period
    # -----------------------------
    if use_angles:
        if physical_to_latent:
            # convert to X, Y, Z
            rho = 1.0

            alpha = pars.pop("alpha", None)
            delta = pars.pop("delta", None)
            period = pars.pop("period", None)

            X = rho * np.cos(delta) * np.cos(alpha)
            Y = rho * np.cos(delta) * np.sin(alpha)
            Z = rho * np.sin(delta)
            pars_out.update({"X": X, "Y": Y, "Z": Z, "period": period})
        else:
            X = pars.pop("X", None)
            Y = pars.pop("Y", None)
            Z = pars.pop("Z", None)
            period = pars.pop("period", None)
            rho = np.sqrt(X**2 + Y**2 + Z**2)
            delta = np.arcsin(Z / rho)
            alpha = np.arctan2(Y, X) % (2 * np.pi)
            pars_out.update({"alpha": alpha, "delta": delta, "period": period})
    else:
        # just forward parameters unchanged
        for k in ["alpha", "delta", "period"]:
            if k in pars:
                pars_out[k] = pars.pop(k)

    # -----------------------------
    # Shape ratios
    # -----------------------------
    if use_shape:
        a_b = pars.pop("a_b", None)
        a_c = pars.pop("a_c", None)
        if physical_to_latent:
            u_a_b = logit((a_b - 1) / 4)
            u_a_c = logit((a_c - a_b) / (5 - a_b))
            pars_out.update({"u_a_b": u_a_b, "u_a_c": u_a_c})
        else:
            u_a_b = pars.pop("u_a_b", None)
            u_a_c = pars.pop("u_a_c", None)
            a_b = 4 * sigmoid(u_a_b) + 1
            a_c = (5 - a_b) * sigmoid(u_a_c) + a_b
            pars_out.update({"a_b": a_b, "a_c": a_c})
    else:
        # just forward parameters unchanged
        for k in ["a_b", "a_c"]:
            if k in pars:
                pars_out[k] = pars.pop(k)

    # -----------------------------
    # Phase
    # -----------------------------
    if use_phase:
        phi0 = pars.pop("W0", None)
        if physical_to_latent:
            u_phi0 = logit((phi0 + np.pi / 2) / np.pi)
            pars_out["u_W0"] = u_phi0
        else:
            u_phi0 = pars.pop("u_W0", None)
            phi0 = np.pi * sigmoid(u_phi0) - np.pi / 2
            pars_out["W0"] = phi0
    else:
        pars_out["W0"] = pars["W0"]

    # -----------------------------
    # Filter-dependent parameters
    # -----------------------------
    if use_filter_dependent:
        H_keys = sorted([k for k in pars if k.startswith("H")])
        if physical_to_latent:
            key_read = ""
            key_write = "u_"
        else:
            key_read = "u_"
            key_write = ""

        for H_key in H_keys:
            band = H_key[1:]

            G1_key_read = f"{key_read}G1{band}"
            G2_key_read = f"{key_read}G2{band}"
            G1_key_write = f"{key_write}G1{band}"
            G2_key_write = f"{key_write}G2{band}"

            H = pars.pop(H_key)
            G1 = pars.pop(G1_key_read)
            G2 = pars.pop(G2_key_read)

            if physical_to_latent:
                u_H = H
                u_G1 = sc_logit(G1)
                L, U = compute_LU_bounds(G1)
                u_G2 = logit((G2 - L) / (U - L))
                pars_out.update({H_key: u_H, G1_key_write: u_G1, G2_key_write: u_G2})
            else:
                u_H = H
                u_G1 = G1
                u_G2 = G2
                H_phys = u_H
                G1_phys = sc_sigmoid(u_G1)
                L, U = compute_LU_bounds(G1_phys)
                G2_phys = L + (U - L) * sigmoid(u_G2)
                pars_out.update(
                    {H_key: H_phys, G1_key_write: G1_phys, G2_key_write: G2_phys}
                )
    else:
        # just keep remaining keys
        pars_out.update(pars)

    # -----------------------------
    # reference epoch
    # -----------------------------
    if "t0" in pars:
        pars_out["t0"] = pars["t0"]

    return pars_out


def lmfit_to_dict(params):
    """
    
    """
    return {name: p.value for name, p in params.items()}


def dict_to_lmfit(par_dict, template_params):
    """
    
    """

    new_params = lmfit.Parameters()

    for name, value in par_dict.items():
        if name in template_params:
            p = template_params[name]
            new_params.add(
                name,
                value=value,
                min=p.min,
                max=p.max,
                vary=p.vary,
            )
        else:
            new_params.add(name, value=value)

    return new_params
