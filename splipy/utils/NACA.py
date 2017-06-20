# -*- coding: utf-8 -*-

from splipy import Curve, BSplineBasis
import splipy.surface_factory as SurfaceFactory
import numpy as np

__all__ = ['camber', 'NACA']


def camber(M, P, order=5):
    """ Create the NACA centerline used for wing profiles. This is given as
    an exact quadratic piecewise polynomial y(x),
    see http://airfoiltools.com/airfoil/naca4digit. The method will produce
    one of two representations: For order<5 it will create x(t)=t and
    for order>4 it will create x(t) as qudratic in t and stretched towards the
    endpointspoints creating a more optimized parametrization.
    :param M: Max camber height (y) given as percentage 0% to 9% of length
    :type  M: Int 0<M<10
    :param P: Max camber position (x) given as percentage 0% to 90% of length
    :type  P: Int 0<P<10
    :return : Exact centerline representation
    :rtype  : Curve
    """
    # parametrized by x=t or x="t^2" if order>4
    M = M / 100.0
    P = P / 10.0
    basis = BSplineBasis(order)
    # basis.insert_knot([P]*(order-2)) # insert a C1-knot
    for i in range(order - 2):
        basis.insert_knot(P)

    t = basis.greville()
    n = len(t)
    x = np.zeros((n, 2))
    for i in range(n):
        if t[i] <= P:
            if order > 4:
                x[i, 0] = t[i]**2 / P
            else:
                x[i, 0] = t[i]
            x[i, 1] = M / P / P * (2 * P * x[i, 0] - x[i, 0] * x[i, 0])
        else:
            if order > 4:
                x[i, 0] = (t[i]**2 - 2 * t[i] + P) / (P - 1)
            else:
                x[i, 0] = t[i]
            x[i, 1] = M / (1 - P) / (1 - P) * (1 - 2 * P + 2 * P * x[i, 0] - x[i, 0] * x[i, 0])
    N = basis.evaluate(t)
    controlpoints = np.linalg.solve(N, x)
    return Curve(basis, controlpoints)


def NACA(M, P, X, n=40, order=5, closed=False):
    """ Create the NACA 4 digit airfoil. This is generated as an approximation
    through the use of SurfaceFactory.thicken functions.
    :param M: Max camber height (y) given as percentage 0% to 9% of length
    :type  M: Int 0<M<10
    :param P: Max camber position (x) given as percentage 0% to 90% of length
    :type  P: Int 0<P<10
    :param X: Thickness given in percentage 1
    :type  X: Int 0<X<40
    :param n: Number of knot spans in resulting parametrization
    :type  n: Int
    :param order: Order of resulting spline curve
    :type  order: Int
    :return : Approximation of NACA wing profile
    :rtype  : Curve
    """

    n = int(n / 4.0)  # knot spans on each side
    center_line = camber(M, P, order)
    new_knots1 = np.linspace(0, P / 10.0, n + 2)
    new_knots2 = np.linspace(P / 10.0, 1, n + 2)
    new_knots = np.zeros(2 * n)
    new_knots[:n] = new_knots1[1:-1]
    new_knots[n:] = new_knots2[1:-1]
    center_line.insert_knot(new_knots)
    T = X / 100.0

    def thickness(x):
        a0 = 0.2969
        a1 = -0.126
        a2 = -0.3516
        a3 = 0.2843
        a4 = -0.1036 if closed else -0.1015
        return T / 0.2 * (a0 * np.sqrt(x) + a1 * x + a2 * x**2 + a3 * x**3 + a4 * x**4)

    surf = SurfaceFactory.thicken(center_line, thickness)
    _, _, top, btm = surf.edges()
    top.reverse()
    top.append(btm)

    if closed:
        top = top.make_periodic(0)

    return top
