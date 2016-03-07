.. Splipy documentation master file, created by
   sphinx-quickstart on Sun Jan 24 12:47:58 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

======
Splipy
======

Splipy is a python package for building spline geometries. It has been designed
primarily for analysis use, and therefore allows fine-grained control over many
aspects which is not possible to achieve with conventional CAD tools.

For a quick start, please see the :doc:`tutorial`.


Overview
========

The *splipy* module exports a handful of convenient classes. Please see
:doc:`basic_classes` for more information.

.. currentmodule:: splipy
.. autosummary::
   :nosignatures:

   BSplineBasis
   SplineObject
   Curve
   Surface
   Volume

Some convenience utilities for commonly seen objects can be found in the factory
modules. Please see :doc:`factories` for more information.

.. autosummary::

   splipy.curve_factory
   splipy.surface_factory
   splipy.volume_factory

Miscellaneous utilities can be found in the ``splipy.utils`` module.

.. autosummary::

   splipy.utils.curve
   splipy.utils.image
   splipy.utils.refinement
   splipy.utils.smooth

For input and output, please refer to :doc:`io`.


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

