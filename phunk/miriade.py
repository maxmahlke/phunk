from astropy.time import Time
import pandas as pd
import requests


def query(name, epochs):
    """Gets asteroid ephemerides from IMCCE Miriade.

    Parameters
    ----------
    name : str
        Name or designation of asteroid.
    epochs : list
        List of observation epochs in MJD format.

    Returns
    -------
    pd.DataFrame - Input dataframe with ephemerides columns appended
                False - If query failed somehow
    """

    # Pass sorted list of epochs to speed up query
    # Have to convert them to JD
    epochs = [Time(str(e), format="mjd").jd for e in epochs]
    files = {"epochs": ("epochs", "\n".join(["%.6f" % epoch for epoch in epochs]))}

    # ------
    # Query Miriade for phase angles
    url = "http://vo.imcce.fr/webservices/miriade/ephemcc_query.php"

    params = {
        "-name": f"a:{name}",
        "-mime": "json",
        "-tcoor": "5",
        "-output": "--jul",
        "-tscale": "UTC",
    }

    # Execute query
    try:
        r = requests.post(url, params=params, files=files, timeout=50)
    except requests.exceptions.ReadTimeout:
        return False
    j = r.json()

    # Read JSON response
    try:
        ephem = pd.DataFrame.from_dict(j["data"])
    except KeyError:
        return False

    return ephem
