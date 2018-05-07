================
Input and output
================

Splipy supports reading from and writing to a number of different file formats.
This is handled through a number of classes in the ``splipy.io`` module, each of
which have identical interfaces (although not all formats support all
functions).

The file objects are intended to be used as context managers, i.e.

.. code:: python

   with G2('some_file.g2') as f:
       f.write(...)
       objs = f.read()

The available formats are:

- :class:`splipy.io.G2`: G2-format (as used by GoTools)
- :class:`splipy.io.SVG`: Scalable Vector Graphics
- :class:`splipy.io.STL`: Stereolitography format
- :class:`splipy.io.OpenFOAM`: OpenFOAM mesh format

The G2, SVG and STL readers can read and write spline objects directly. The OpenFOAM
file format only writes :class:`splipy.SplineModel.SplineModel` objects.


Common interface
================

.. autoclass:: splipy.io.master.MasterIO
   :members:
   :special-members: __init__
