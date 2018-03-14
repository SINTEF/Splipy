from splipy import Curve, Surface, Volume, BSplineBasis
import numpy as np
import unittest
import pytest

def get_spline(spline, n, p, rational=False):
    basis = BSplineBasis(p, [0]*(p-1) + list(range(n-p+2)) + [n-p+1]*(p-1))
    if spline is 'curve':
        return Curve(basis, rational=rational)
    elif spline is 'surface':
        return Surface(basis, basis, rational=rational)
    elif spline is 'volume':
        return Volume(basis, basis, basis, rational=rational)
    return None

def knot_insert(what, spline_n, spline_p, insert_n):
    spline = get_spline(what, spline_n, spline_p)
    new_knots = np.linspace(spline.start(0), spline.end(0), insert_n, endpoint=False)
    spline.insert_knot(new_knots, 0)

def raise_order(what, spline_n, spline_p, raise_n):
    spline = get_spline(what, spline_n, spline_p)
    amount = [raise_n]*spline.pardim
    spline.raise_order(*amount)

def move_inplace(spline):
    spline += (1,1,1)

def move_function(spline):
    spline.translate((1,1,1))

def move_operator(spline):
    spline = spline + (1,1,1)

def grid_evaluate(spline, n):
    # create n x...x n evaluation grid
    u = tuple([np.linspace(u0, u1, n) for u0,u1 in zip(spline.start(), spline.end())])
    x = spline(*u) # u has 1, 2 or 3 components

def derivative_evaluate(spline, n):
    u = tuple([np.linspace(u0, u1, n) for u0,u1 in zip(spline.start(), spline.end())])
    d = [0]*spline.pardim
    d[0] = 1
    x = spline.derivative(*u,d=d) # u has 1, 2 or 3 components

def tangent_evaluate(spline, n):
    u = tuple([np.linspace(u0, u1, n) for u0,u1 in zip(spline.start(), spline.end())])
    x = spline.tangent(*u) # u has 1, 2 or 3 components




@pytest.mark.benchmark(group="knot-insert")
def test_curve_ki(benchmark):
    benchmark(knot_insert, 'curve', 350, 3, insert_n=5)

@pytest.mark.benchmark(group="knot-insert")
def test_surface_ki(benchmark):
    benchmark(knot_insert, 'surface', 110, 3, insert_n=5)

@pytest.mark.benchmark(group="knot-insert")
def test_volume_ki(benchmark):
    benchmark(knot_insert, 'volume', 40, 3, insert_n=5)





@pytest.mark.benchmark(group="raise-order")
def test_curve_rs(benchmark):
    benchmark(raise_order, 'curve', 150, 3, raise_n=1)

@pytest.mark.benchmark(group="raise-order")
def test_surface_rs(benchmark):
    benchmark(raise_order, 'surface', 80, 3, raise_n=1)

@pytest.mark.benchmark(group="raise-order")
def test_volume_rs(benchmark):
    benchmark(raise_order, 'volume', 30, 3, raise_n=1)




@pytest.mark.benchmark(group="move")
def test_move_inplace(benchmark):
    spline = get_spline('volume', 30, 3)
    benchmark(move_inplace, spline)

@pytest.mark.benchmark(group="move")
def test_move_function(benchmark):
    spline = get_spline('volume', 30, 3)
    benchmark(move_function, spline)

@pytest.mark.benchmark(group="move")
def test_move_operator(benchmark):
    spline = get_spline('volume', 30, 3)
    benchmark(move_operator, spline)




@pytest.mark.benchmark(group="evaluate-splineobject")
def test_eval(benchmark):
    spline = get_spline('volume', 30, 3)
    benchmark(grid_evaluate, spline, 15)

@pytest.mark.benchmark(group="evaluate-splineobject")
def test_eval_rational(benchmark):
    spline = get_spline('volume', 30, 3, True)
    benchmark(grid_evaluate, spline, 15)

@pytest.mark.benchmark(group="evaluate-splineobject")
def test_deriv(benchmark):
    spline = get_spline('volume', 30, 3)
    benchmark(derivative_evaluate, spline, 15)

@pytest.mark.benchmark(group="evaluate-splineobject")
def test_deriv_rational(benchmark):
    spline = get_spline('volume', 30, 3, True)
    benchmark(derivative_evaluate, spline, 15)

@pytest.mark.benchmark(group="evaluate-splineobject")
def test_tangent(benchmark):
    spline = get_spline('volume', 30, 3)
    benchmark(tangent_evaluate, spline, 15)

@pytest.mark.benchmark(group="evaluate-splineobject")
def test_tangent_rational(benchmark):
    spline = get_spline('volume', 30, 3, True)
    benchmark(tangent_evaluate, spline, 15)
