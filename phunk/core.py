import numpy as np

import pandas as pd
import rocks

import phunk


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

        self.phase = np.array(phase)
        self.mag = np.array(mag)
        self.mag_err = np.array(mag_err)
        self.target = target
        self.epoch = np.array(epoch)
        self.ra = np.array(ra)
        self.dec = np.array(dec)
        self.band = np.array(band)

        if target is not None:
            self.target = rocks.Rock(target)

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

    def fit(self, models=None, weights=None):
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
                getattr(self, model).fit(self.phase, self.mag, weights)

        # if "HG1G2" in models:
        #     # HG1G2
        #     for band in self.bands:
        #         hg1g2 = HG1G2(band=band)
        #         hg1g2.fit(
        #             phase=getattr(self, f"obs_{band}").phase.values,
        #             mag=getattr(self, f"obs_{band}").mred.values,
        #             mag_err=getattr(self, f"obs_{band}").dm.values,
        #         )
        #         setattr(self, f"HG1G2_{band}", hg1g2)
        #
        # if "sHG1G2" in models:
        #     # sHG1G2
        #     shg1g2 = sHG1G2(obs=self.obs)
        #     shg1g2.fit_mm(
        #         phase=self.obs.phase.values,
        #         mag=self.obs.mred.values,
        #         mag_err=self.obs.dm.values,
        #         # mred=self.obs.mred.values,
        #         # mred_err=self.obs.dm.values,
        #         ra=self.obs.ra.values,
        #         dec=self.obs.dec.values,
        #         # filters=self.obs.band.values,
        #         bands=self.obs.band.values,
        #     )
        #     self.sHG1G2 = shg1g2

    # def save(self, db=None, models=None):
    #     """Save fit parameters in atlas phase curve database."""
    #
    #     if db is None:
    #         data = atlas.data.load_fits()
    #     elif isinstance(db, atlas.data.Database):
    #         data = db.data
    #     else:
    #         raise ValueError("Invalid db type.")
    #
    #     if models is None:
    #         models = ["HG1G2", "sHG1G2"]
    #
    #     entry = {}
    #
    #     # Store target parameters
    #     entry["name"] = self.target.name
    #     # entry["number"] = self.target.number
    #     # entry["class_"] = self.target.class_
    #     # entry["D"] = self.target.diameter.value
    #     # entry["taxonomy"] = self.target.taxonomy.class_.value
    #     # entry["complex"] = self.target.taxonomy.complex.value
    #     # entry["pV"] = self.target.albedo.value
    #
    #     # Store phase curve parameters
    #     # entry["N"] = len(self.obs)
    #
    #     for band in ["o", "c"]:
    #         entry[f"N_{band}_fit"] = getattr(self, f"N_{band}")
    #         entry[f"phase_min_{band}_fit"] = getattr(self, f"phase_min_{band}")
    #         entry[f"phase_max_{band}_fit"] = getattr(self, f"phase_max_{band}")
    #
    #     # Store model parameters
    #     if "HG1G2" in models:
    #         for band in self.bands:
    #             for param in ["H", "G1", "G2"]:
    #                 entry[f"HG1G2_{band}_{param}"] = getattr(
    #                     getattr(self, f"HG1G2_{band}"), param
    #                 )
    #                 entry[f"HG1G2_{band}_{param}_err"] = getattr(
    #                     getattr(self, f"HG1G2_{band}"), f"{param}_err"
    #                 )
    #     if "sHG1G2" in models:
    #         for band in self.bands:
    #             for param in ["H", "G1", "G2"]:
    #                 entry[f"sHG1G2_{band}_{param}"] = getattr(
    #                     getattr(self, "sHG1G2"), f"{param}_{band}"
    #                 )
    #                 entry[f"sHG1G2_{band}_{param}_err"] = getattr(
    #                     getattr(self, "sHG1G2"), f"{param}_{band}_err"
    #                 )
    #
    #         entry[f"sHG1G2_R"] = getattr(self, "sHG1G2").R
    #         entry[f"sHG1G2_alpha0"] = getattr(self, "sHG1G2").alpha0
    #         entry[f"sHG1G2_delta0"] = getattr(self, "sHG1G2").delta0
    #         entry[f"sHG1G2_R_err"] = getattr(self, "sHG1G2").R_err
    #         # entry[f"sHG1G2_alpha0_err"] = getattr(self, "sHG1G2").delta0_err
    #         # entry[f"sHG1G2_delta0_err"] = getattr(self, "sHG1G2").alpha0_err
    #
    #     # Concatenate or overwrite in database
    #     if entry["name"] in data.name.values:
    #         data = data[data.name != entry["name"]]
    #     entry["is_fit"] = True
    #     # data = data.append(entry, ignore_index=True)
    #     data = pd.concat([data, pd.DataFrame(data=entry, index=[0])], ignore_index=True)
    #
    #     if db is None:
    #         atlas.data.save_fits(data)
    #     else:
    #         db.data = data
