import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

import phunk
import rocks


def get_colors(N, cmap="turbo"):
    """
    Get a list of unique colors.

    Parameters
    ----------
    N : int
        The number of unique colors to return.
    cmap : str
        The matplotlib colormap to sample. Default is 'turbo'

    Returns
    -------
    list of str
        A list of color-hexcodes.

    """
    COLORS = plt.get_cmap(cmap, N)
    return [mpl.colors.rgb2hex(COLORS(i)[:3]) for i in range(N)]


def plot_pc(pc, models, band=None, save=None):
    """Plot phase curve and model fits."""

    fig, ax = plt.subplots()

    if band is None:
        band = pc.band[0]

    phase = pc.phase[pc.band == band]
    mag = pc.mag[pc.band == band]
    mag_err = pc.mag_err[pc.band == band]

    # ------
    # Plot observations
    ax.errorbar(
        phase,
        mag,
        yerr=mag_err,
        ls="",
        marker="x",
        c="black",
        label=f"Observations in {band}",
    )

    # ------
    # Plot model fits
    COLORS = get_colors(len(models))

    for i, model in enumerate(models):
        model = getattr(pc, model)  # switch to actual model instance

        phase_eval = np.linspace(0, pc.phase.max(), 100)
        params = model.PARAMS

        # You have to specify the band if plotting sHG1G2
        if model.NAME == "sHG1G2":
            mag_eval = phunk.models.HG1G2().eval(
                phase_eval,
                H=getattr(model, f"H_{band}"),
                G1=getattr(model, f"G1_{band}"),
                G2=getattr(model, f"G2_{band}"),
            )
            params = [f"H_{band}", f"G1_{band}", f"G2_{band}"]  # + params

        else:
            mag_eval = model.eval(phase_eval)

        label = ": ".join(
            [
                model.NAME,
                ", ".join(
                    [
                        f"{p}: {v:.2f}"
                        for p, v in zip(params, [getattr(model, p) for p in params])
                    ]
                ),
            ]
        )

        ax.plot(phase_eval, mag_eval, c=COLORS[i], label=label, ls="-")

    # ------
    # Set up figure
    ax.yaxis.set_inverted(True)
    ax.set(
        xlabel="Phase Angle / deg", ylabel="Reduced Magnitude", xlim=(0, pc.phase.max())
    )

    ax.legend(
        title=f"({pc.target.number}) {pc.target.name}"
        if isinstance(pc.target, rocks.Rock)
        else None
    )

    if save is None:
        plt.show()
    else:
        fig.savefig(save)
        print(f"Saved figure under {save}")
