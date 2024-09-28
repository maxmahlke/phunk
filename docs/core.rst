Basic Usage
-----------

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

   Where applicable, ``phunk`` expects values to be in degrees rather than radians.

Providing ``phase`` and ``mag`` is enough to fit most :ref:`models <models>`.

Multi-band Phase Curves
+++++++++++++++++++++++

You can pass observations in different bands to the ``PhaseCurve`` class by making use of the ``band`` argument,
which should be a list which defines the band for each observation.

.. code-block:: python

   phase = [0.57, 1.09, 3.20, 10.99, 14.69, 20.42]
   mag = [6.555, 6.646, 6.793, 7.130, 7.210, 7.414]
   band = ["V", "V", "G", "G", "R", "V"]

   pc = PhaseCurve(phase=phase, mag=mag, band=band)

At any point, you can check the bands of your observations using the ``bands`` attribute of the phase curve.

.. code-block:: python

   >>> pc.bands
   ["G", "R", "V"]

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

Datapoints can be weighted by providing the ``weights`` argument.

.. code-block:: python

   pc.fit("HG1G2", weights=[...])

Once the models have been fit, you can access the model parameters as attributes of
the ``PhasCurve`` via the dot notation.

.. code-block:: python

   pc.HG.H
   pc.HG1G2.G1
   pc.HG1G2.G1_err
   pc.sHG1G2.R

All available model attributes are given in the model description.

Fitting multi-band Curves
+++++++++++++++++++++++++

When fitting multi-band observations, ``phunk`` automatically separates the observations by band and fits the separate
phase curves. The model parameters of the different fitted phase curves get the band name as suffix.

.. code-block:: python

   epoch = [35193, 35194, 35198, 35214, 35223, 35242]
   mag = [6.555, 6.646, 6.793, 7.130, 7.210, 7.414]
   band = ["o", "o", "o", "c", "c", "c"]

   pc = PhaseCurve(mag=mag, epoch=epoch, band="o")
   pc.fit(['HG1G2'])  # fits two phase curves, one in "c" and the other in "o"

   pc.HG1G2.Hc  # absolute magnitude in "c"
   pc.HG1G2.Ho  # absolute magnitude in "o"
   pc.HG1G2.H   # undefined

.. Note::

   The asteroid-specific parameters ``alpha``, ``delta``, and ``R`` of the ``sHG1G2`` remain without band-specific suffix
   as they are uniquely defined, independent of observation band.

Constrained solutions
+++++++++++++++++++++

The ``G1`` and ``G2`` parameters of the ``HG1G2`` and ``sHG1G2`` :ref:`models <models>` can be constrained to yield
only physically-meaningful solutions (decreasing brightness with increasing phase angle). By default, they are limited
to values between ``0`` and ``1``. If you set ``constrain_g1g2=True`` when fitting ``HG1G2`` and ``sHG1G2``, the ``1-G1-G2>=0``
inequality constraint is further applied.

.. code-block:: python

   pc.fit("sHG1G2", constrain_g1g2=True)

Note that this might not be advantageous
for phase curve fits: The solution will still be bad, and if unconstrained, you can use ``G1`` and ``G2`` to quickly identify failed fits.



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
