Basic Usage
-----------

.. dot graph
..
.. phase curve -> models -> a + b exp c
.. phase curve -> models -> HG1G2
.. phase curve -> models -> HG12
.. phase curve -> models -> HG
.. phase curve -> models -> HG12*
.. phase curve -> models -> sHG1G2
..
.. phase curve -> phase
.. phase curve -> mag
.. phase curve -> band
.. phase curve -> redmag
.. phase curve -> epoch
..
.. phase curve -> target -> rocks.Rock
.. phase curve -> target -> ephemerides
.. phase curve -> target -> ephemerides -> phase
.. phase curve -> target -> ephemerides -> ra
.. phase curve -> target -> ephemerides -> dec
..
.. Specify epoch in MJD
..
.. pc = phunk.PhaseCurve(mag=..., target='Massalia', epoch=...)
.. pc.fit(sHG1G2)
.. pc.plot()

Providing Observations
======================

The first step is to create a ``PhaseCurve`` instance, representing the observed phase curve
of the target. In its simplest form, this consists of providing a list of ``phase`` (in degrees)
and a list of ``mag`` to the ``PhaseCurve``.

.. code-block:: python

   from phunk import PhaseCurve

   phase = [0.57, 1.09, 3.20, 10.99, 14.69, 20.42]
   mag = [6.555, 6.646, 6.793, 7.130, 7.210, 7.414]

   pc = PhaseCurve(phase=phase, mag=mag)

.. Note::

  ``phunk`` uses ``mag`` to refer to the **reduced** magnitude.

.. Note::

   If applicable, ``phunk`` expects values to be in degrees rather than radians.

Providing ``phase`` and ``mag`` is enough to fit most :ref:`models <models>`.

Computing Ephemerides
+++++++++++++++++++++

If you do not have the ``phase`` angle or you require your target's ephemerides to fit the :ref:`sHG1G2 <sHG1G2>` model, ``phunk``
can query the required values for you from the IMCCE's `Miriade <https://ssp.imcce.fr/webservices/miriade/>`_ webservice.
For this, you have to specify the ``epoch`` of observation in MJD format and an identifier of your target. ``phunk`` uses `rocks <https://rocks.readthedocs.io>`_
to resolve target identities.

.. code-block:: python

   from phunk import PhaseCurve

   epoch = [35193, 35194, 35198, 35214, 35223, 35242]
   mag = [6.555, 6.646, 6.793, 7.130, 7.210, 7.414]

   pc = PhaseCurve(mag=mag, epoch=epoch, target="massalia")

Fitting Models
==============

To fit one of the available :ref:`photometric models <models>`, use the ``.fit`` method of the ``PhaseCurve``
and provide a list of models to fit.

.. code-block:: python

   pc.fit(["HG1G2"])
   pc.fit(["sHG1HG", "HG1G2"])

If you don't provide any argument, ``phunk`` will fit all implemented models.

.. Datapoints can be weighted by providing the ``weights`` argument.
..
.. .. code-block:: python
..
..    pc.fit("HG1G2", weights=[...])

Accessing results
=================

Once the models have been fit, you can access the model parameters as attributes of
the ``PhasCurve`` via the dot notation.

.. code-block:: python

   pc.HG.H
   pc.HG1G2.G1
   pc.HG1G2.G1_err
   pc.sHG1G2.R

All available model attributes are given in the model description.

Plotting results
================

Use the ``.plot`` method of the ``PhaseCurve`` class to plot phase curves.
You can select which models to add to the plot using the ``models`` argument.
The plot will open in an interactive window by default. Provide a path to the ``save``
argument to save the plot under the specified path.

.. code-block:: python

   pc.plot()
   pc.plot(models=["sHG1G2"])
   pc.plot(models=["sHG1G2"], save="graphics/massalia_gehrels_shg1g2.png")

.. Note::

   You need to ``fit`` a model before you can ``plot`` it.

Multi-band phase curves
=======================
pc.bands
