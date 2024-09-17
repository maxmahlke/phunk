import phunk
import numpy as np
import rocks


class PhaseCurve:
    """Phase curve of a given asteroid."""

    def __init__(
        self,
        phase=None,
        mag=None,
        mag_err=None,
        target=None,
        epoch=None,
        ra=None,
        dec=None,
        band=None,
    ):
        """Create a phase curve and add observations to it.

        Parameters
        ----------
        target : rocks.Rock
            The target asteroid.
        obs : pd.DataFrame
            Observations describing the phase curve. Relevant columns are

            phase (degrees), mred, mred_err, band
        """

        # Data
        self.phase = np.array(phase)
        self.mag = np.array(mag)

        self.mag_err = (
            np.array(mag_err) if mag_err is not None else np.zeros(self.mag.shape)
        )

        # Metadata
        self.target = target

        self.epoch = np.array(epoch) if epoch is not None else epoch
        self.ra = np.array(ra) if ra is not None else ra
        self.dec = np.array(dec) if dec is not None else dec
        self.band = np.array(band).astype(str) if band is not None else band

        if self.band is None:
            print("No observation bands provided. Assuming they are all 'V' band.")
            self.band = np.array(["V"] * len(self.phase))

        if target is not None:
            self.target = rocks.Rock(target)

        self.fitted_models = set()  # keep track of models fit to data

        # TODO: Make these properties @property methods
        #
        # self.bands = obs.band.unique()
        #
        # for band in self.bands:
        #     obs_ = obs[obs.band == band]
        #     setattr(self, f"obs_{band}", obs_)
        #     setattr(self, f"N_{band}", len(obs_))
        #     setattr(self, f"phase_min_{band}", obs_.phase.min())
        #     setattr(self, f"phase_max_{band}", obs_.phase.max())

    def fit(self, models=None):
        """Fit the phase curve in the different bands with the different models."""

        if models is None:
            models = phunk.models.MODELS

        for model in models:
            if model not in phunk.models.MODELS:
                raise ValueError(
                    f"Unknown model '{model}'. Expected one of {phunk.models.MODELS}"
                )

            setattr(self, model, getattr(phunk.models, model)())

            if getattr(self, model).is_fittable(self):
                getattr(self, model).fit(self)

    def get_ephems(self):
        """Query ephemerides of target at time of observation.

        Note
        ----
        Sets the 'phase', 'ra', and 'dec' attributes. Requires internet connection.
        """
        print("Querying ephemerides via IMCCE Miriade..")
        ephem = phunk.miriade.query(self.target.name, self.epoch)

        self.phase = ephem["phase"]
        self.ra = np.degrees(ephem["ra_j2000"])
        self.dec = np.degrees(ephem["dec_j2000"])

    def plot(self, models=None, band=None, save=None):
        """Plot phase curve and model fits.

        Parameters
        ----------
        models : list of str
            Name of models to plot. By default, all fitted models are plotted.
        band : str
            If plotting sHG1G2 model, observation band to plot
        """
        if models is None:
            models = self.fitted_models

        phunk.plotting.plot_pc(self, models, band, save)
