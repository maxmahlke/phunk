Basic Usage
-----------

dot graph

phase curve -> models -> a + b exp c
phase curve -> models -> HG1G2
phase curve -> models -> HG12
phase curve -> models -> HG
phase curve -> models -> HG12*
phase curve -> models -> sHG1G2

phase curve -> phase
phase curve -> mag
phase curve -> band
phase curve -> redmag
phase curve -> epoch

phase curve -> target -> rocks.Rock
phase curve -> target -> ephemerides
phase curve -> target -> ephemerides -> phase
phase curve -> target -> ephemerides -> ra
phase curve -> target -> ephemerides -> dec

Providing Observations
======================

The first step is to create a ``PhaseCurve`` instance, representing the observed phase curve
of the target. In its simplest form, this consists of providing a list of ``phase`` (in degrees)
and a list of ``mag`` to the ``PhaseCurve``.

.. code-block:: python

   from phunk import PhaseCurve

   phase = []
   mag = []

   pc = PhaseCurve(phase=phase, mag=mag)

.. Note::

  mag
  mag_app
  mag_red

Computing Ephemerides
+++++++++++++++++++++

Specify the target

Fitting Models
==============

.. TODO: Link to models page

To fit one of the available photometric models, use the ``.fit`` method of the ``PhaseCurve``
and provide the name of the model or a list of models to fit.

.. code-block:: python

   pc.fit("HG1G2")
   pc.fit(["sHG1HG", HG1G2"])

If you don't provide any argument, ``phunk`` will fit all models.
Datapoints can be weighted by providing the ``weights`` argument.

.. code-block:: python

   pc.fit("HG1G2", weights=[...])

Accessing results
=================

Once the models have been fit, you can access the model parameters as attributes of
the ``PhasCurve`` via the dot notation.

.. code-block:: python

   pc.HG1G2.H
   pc.HG1G2.G1
   pc.HG1G2.G2

All available model attributes are given in the model description.

Plotting results
================
