=============
Basic classes
=============


BSplineBasis
============

.. autoclass:: splipy.BSplineBasis
   :members:
   :special-members: __init__, __len__, __getitem__
   :inherited-members:


SplineObject
============

.. autoclass:: splipy.SplineObject
   :members:
   :special-members: __init__


Curve
=====

.. autoclass:: splipy.Curve
   :members:
   :special-members: __init__, __len__, __getitem__, __setitem__
   :exclude-members: evaluate
   :inherited-members:
   :show-inheritance:

   .. method:: evaluate(u)

      Evaluate the curve at the given parametric values.

      This function returns an *n* × *dim* array, where *n* is the number of
      evaluation points, and *dim* is the physical dimension of the curve.

      If there is only one evaluation point, a vector of length *dim* is
      returned instead.

      :param u: Parametric coordinates in which to evaluate
      :type u: float or [float]
      :return: Geometry coordinates
      :rtype: numpy.array


Surface
=======

.. autoclass:: splipy.Surface
   :members:
   :special-members: __init__, __len__, __getitem__, __setitem__
   :exclude-members: evaluate, evaluate_derivative
   :inherited-members:
   :show-inheritance:

   .. method:: evaluate(u, v)

      Evaluate the surface at the given parametric values.

      This function returns an *n1* × *n2* × *dim* array, where *ni* is the
      number of evaluation points in each direction, and *dim* is the dimension
      of the surface.

      If there is only one evaluation point, a vector of length *dim* is
      returned instead.

      :param u: Parametric coordinates in the first direction
      :type u: float or [float]
      :param u: Parametric coordinates in the second direction
      :type u: float or [float]
      :return: Geometry coordinates
      :rtype: numpy.array

   .. method:: evalute_derivative(u, v, [d=(1,1)])

      Evaluate the derivative of the surface at the given parametric values.

      This function returns an *n1* × *n2* × *dim* array, where *ni* is the
      number of evaluation points in direction *i*, and *dim* is the dimension
      of the surface.

      If there is only one evaluation point, a vector of length *dim* is
      returned instead.

      :param u: Parametric coordinates in the first direction
      :type u: float or [float]
      :param u: Parametric coordinates in the second direction
      :type u: float or [float]
      :param (int) d: Order of derivative to compute
      :return: Derivatives
      :rtype: numpy.array


Volume
======

.. autoclass:: splipy.Volume
   :members:
   :special-members: __init__, __len__, __getitem__, __setitem__
   :exclude-members: evaluate, evaluate_derivative
   :inherited-members:
   :show-inheritance:

   .. method:: evaluate(u, v, w)

      Evaluate the volume at the given parametric values.

      This function returns an *n1* × *n2* × *n3* × *dim* array, where *ni* is
      the number of evaluation points in each direction, and *dim* is the
      dimension of the volume.

      If there is only one evaluation point, a vector of length *dim* is
      returned instead.

      :param u: Parametric coordinates in the first direction
      :type u: float or [float]
      :param v: Parametric coordinates in the second direction
      :type v: float or [float]
      :param w: Parametric coordinates in the third direction
      :type w: float or [float]
      :return: Geometry coordinates
      :rtype: numpy.array

   .. method:: evalute_derivative(u, v, w, [d=(1,1,1)])

      Evaluate the derivative of the volume at the given parametric values.

      This function returns an *n1* × *n2* × *n3* × *dim* array, where *ni* is
      the number of evaluation points in direction *i*, and *dim* is the
      dimension of the volume.

      If there is only one evaluation point, a vector of length *dim* is
      returned instead.

      :param u: Parametric coordinates in the first direction
      :type u: float or [float]
      :param v: Parametric coordinates in the second direction
      :type v: float or [float]
      :param w: Parametric coordinates in the third direction
      :type w: float or [float]
      :param (int) d: Order of derivative to compute
      :return: Derivatives
      :rtype: numpy.array
