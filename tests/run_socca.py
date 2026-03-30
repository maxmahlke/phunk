import phunk
import pandas as pd

# data
data = pd.read_csv("data/Eugenia_ztf.csv")
cols = ["Date", "cmagpsf"]
sol = pd.read_csv("data/Eugenia_socca.csv")

# Phase
pc = phunk.PhaseCurve(
    target="Eugenia",
    epoch=data["Date"].values - 2400000.5,
    phase=data["Phase"].values,
    mag=data["cmagpsf"].values,
    band=data["cfid"].values,
)


p0 = {
    "H": 8.07602859,
    "G1": 0.53,
    "G2": 0.1746,
    "period": 0.23746434,
    "alpha": 119.987,
    "delta": -24.948,
    "a_b": 1.3422,
    "a_c": 2.2407,
    "W0": -0.57,
    "t0": 2459580.5,
}
pc.fit(models=["SOCCA"], p0=p0)


print(pc.SOCCA.p0)
print(pc.SOCCA.H1)
