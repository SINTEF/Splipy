from splipy import BSplineBasis
import numpy as np
import unittest
import pytest


def evaluate_old(basis, nviz):
    t = np.linspace(0, basis.end(), nviz)
    basis.evaluate_old(t)

def evaluate_cython(basis, nviz):
    t = np.linspace(0, basis.end(), nviz)
    basis.evaluate(t)

def evaluate_dense(basis, nviz):
    t = np.linspace(0, basis.end(), nviz)
    basis.evaluate_old(t, sparse=False)

def evaluate_sparse(basis, nviz):
    t = np.linspace(0, basis.end(), nviz)
    basis.evaluate_old(t, sparse=True)

n = 100
p = 4

@pytest.mark.benchmark(group="evaluate-basis")
def test_evaluate_old(benchmark):
    basis = BSplineBasis(p, [0]*(p-1) + list(range(n-p+2)) + [n-p+1]*(p-1))
    benchmark(evaluate_old, basis, 1000)

@pytest.mark.benchmark(group="evaluate-basis")
def test_evaluate_cython(benchmark):
    basis = BSplineBasis(p, [0]*(p-1) + list(range(n-p+2)) + [n-p+1]*(p-1))
    benchmark(evaluate_cython, basis, 1000)

@pytest.mark.benchmark(group="evaluate-basis")
def test_evaluate_dense(benchmark):
    basis = BSplineBasis(p, [0]*(p-1) + list(range(n-p+2)) + [n-p+1]*(p-1))
    benchmark(evaluate_dense, basis, 1000)

@pytest.mark.benchmark(group="evaluate-basis")
def test_evaluate_sparse(benchmark):
    basis = BSplineBasis(p, [0]*(p-1) + list(range(n-p+2)) + [n-p+1]*(p-1))
    benchmark(evaluate_sparse, basis, 1000)

