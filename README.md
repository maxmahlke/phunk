<p align="center">
  <img width="260" src="https://raw.githubusercontent.com/maxmahlke/phunk/main/docs/gfx/logo_phunk.svg">
</p>

<p align="center">
  <a href="https://github.com/maxmahlke/phunk#features"> Features </a> - <a href="https://github.com/maxmahlke/phunk#install"> Install </a> - <a href="https://github.com/maxmahlke/phunk#documentation"> Documentation </a>
</p>

<div align="center">
  <a href="https://img.shields.io/pypi/pyversions/space-phunk">
    <img src="https://img.shields.io/pypi/pyversions/space-phunk"/>
  </a>
  <a href="https://img.shields.io/pypi/v/space-phunk">
    <img src="https://img.shields.io/pypi/v/space-phunk"/>
  </a>
  <a href="https://readthedocs.org/projects/phunk/badge/?version=latest">
    <img src="https://readthedocs.org/projects/phunk/badge/?version=latest"/>
  </a>
</div>


## Features

Observe the phase curve of an asteroid, ...

``` python
>>> from phunk import PhaseCurve
>>> # Observations of (20) Massalia from Gehrels 1956
>>> phase = [0.57, 1.09, 3.20, 10.99, 14.69, 20.42]
>>> mag = [6.555, 6.646, 6.793, 7.130, 7.210, 7.414]
>>> pc = PhaseCurve(phase=phase, mag=mag)
```

..., fit it in one of multiple photometric models, ....

``` python
>>> pc.fit(["HG", "HG1G2", "sHG1G2"])
```

..., and plot / process the results.

``` python
>>> pc.HG1G2.H
>>> pc.HG12.H
>>> pc.plot()
```

Provide a target to ``PhaseCurve`` to have ``phunk`` compute the required ephemerides for you.

``` python
>>> epoch = [35193, 35194, 35198, 35214, 35223, 35242]  # in MJD
>>> pc = PhaseCurve(epoch=epoch, mag=mag, target='massalia')
>>> pc.fit(['sHG1G2'])  # phunk computes required RA, Dec, and phase at epoch of observation
>>> pc.sHG1G2.H
```

## Install

Install from PyPi using `pip`:

     $ pip install space-phunk

The minimum required `python` version is 3.8.


## Documentation

Check out the documentation at [phunk.readthedocs.io](https://phunk.readthedocs.io/en/latest/).

## Acknowledgements

This package uses the photometric model implementations provided by [sbpy](https://sbpy.readthedocs.io/en/stable).
