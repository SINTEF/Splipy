from bisect import bisect_right, bisect_left
import numpy as np
cimport numpy as np
import copy
cimport cython

# cdef inline int int_max(int a, int b): return a if a >= b else b
# cdef inline unsigned int uint_min(unsigned int a, unsigned int b): return a if a <= b else b
def my_bisect_left(array, value, n):
    cdef unsigned int n0 = 0
    cdef unsigned int i
    while n0+1 < n:
        i=np.floor((n0+n-1)/2)
        if array[i] < value:
            n0 = i+1
        else:
            n = i+1
    return n0

def my_bisect_right(array, value, n):
    cdef unsigned int n0 = 0
    cdef unsigned int i
    while n0+1 < n:
        i=np.floor((n0+n-1)/2)
        if array[i] > value:
            n = i+1
        else:
            n0 = i+1
    return n0

@cython.boundscheck(False) # turn off bounds-checking for entire function
def evaluate(knots_in, order_in, eval_t_in, periodic_in, tolerance_in, derivatives_in=0, from_right=True):
    """  Evaluate all basis functions in a given set of points.

    :param knots_in:       Knot vector
    :param order_in:       Parametrization order (polyonomial degree + 1)
    :param eval_t_in:      The parametric coordinate(s) in which to evaluate
    :param periodic_in:    Periodicity of basis
    :param tolerance_in:   Knot tolerance for detecting end-point-evaluations
    :param derivatives_in: Number of derivatives to compute
    :param from_right:     True if evaluation should be done in the limit from above
    :return: Two tuples of all arguments to the scipy sparse csr_matrix
    """

    # wrap everything into c-type datastructures for optimized performence
    cdef np.ndarray knots   = knots_in
    cdef int periodic       = periodic_in
    cdef unsigned int p     = order_in  # knot vector order
    cdef unsigned int n_all = len(knots) - p  # number of basis functions (without periodicity)
    cdef unsigned int n     = len(knots) - p - (periodic+1)  # number of basis functions (with periodicity)
    cdef unsigned int m     = len(eval_t_in)
    cdef unsigned int d     = derivatives_in
    cdef double start       = knots[p-1]
    cdef double end         = knots[n_all]
    cdef double evalT       ;
    cdef double tol         = tolerance_in
    cdef unsigned int mu    = 0
    cdef double[:] t        = copy.deepcopy(eval_t_in)
    cdef double[:] data     = np.zeros(m*p)
    cdef int[:]    indices  = np.zeros(m*p, dtype='int32')
    cdef int[:]    indptr   = np.arange(0,m*p+1,p, dtype='int32')
    cdef double[:] M        = np.zeros(p)  # temp storage to keep all the function evaluations
    cdef unsigned int k,q,j,i
    cdef bint right

    if periodic >= 0:
        # Wrap periodic evaluation into domain
        for i in range(len(t)):
            if t[i] < start or t[i] > end:
                t[i] = (t[i] - start) % (end - start) + start
    for i in range(len(t)):
        right = from_right
        evalT = t[i]
        # Special-case the endpoint, so the user doesn't need to
        if abs(t[i] - end) < tol:
            right = False
        # Skip non-periodic evaluation points outside the domain
        if t[i] < start or t[i] > end:
            continue

        # mu = index of last non-zero basis function
        if right:
            # mu = my_bisect_right(knots, evalT, n_all+p)
            mu = bisect_right(knots, evalT)
        else:
            # mu = my_bisect_left(knots, evalT, n_all+p)
            mu = bisect_left(knots, evalT)
        mu = min(mu, n_all)

        for k in range(p-1):
            M[k] = 0
        M[p-1] = 1  # the last entry is a dummy-zero which is never used
        for q in range(1, p-d):
            j = p-q-1
            k = mu - q -1
            M[j] = M[j] + M[j + 1] * <double>(knots[k + q + 1] - evalT) / (knots[k + q + 1] - knots[k + 1])
            for j in range(p - q , p-1):
                k = mu - p + j  # 'i'-index in global knot vector (ref Hughes book pg.21)
                M[j] = M[j] * <double>(evalT - knots[k]) / (knots[k + q] - knots[k])
                M[j] = M[j] + M[j + 1] * <double>(knots[k + q + 1] - evalT) / (knots[k + q + 1] - knots[k + 1])
            j = p  - 1
            k = mu - 1
            M[j] = M[j] * <double>(evalT - knots[k]) / (knots[k + q] - knots[k])

        for q in range(p-d, p):
            for j in range(p - q - 1, p):
                k = mu - p + j  # 'i'-index in global knot vector (ref Hughes book pg.21)
                if j != p-q-1:
                    M[j] = M[j] * <double>(q) / (knots[k + q] - knots[k])
                if j != p-1:
                    M[j] = M[j] - M[j + 1] * <double>(q) / (knots[k + q + 1] - knots[k + 1])


        for j,k in enumerate(range(i*p, (i+1)*p)):
            data[k]    = M[j]
            indices[k] = (mu-p+j) % n
    return (data, indices, indptr), (m,n)


@cython.boundscheck(False) # turn off bounds-checking for entire function
def snap(knots_in, eval_t_in, tolerance_in):
    """  Snap evaluation points to knots if they are sufficiently close
    as given in by state.state.knot_tolerance. This will modify the input vector t
    :param knots:        Knot vector
    :param eval_t_in:    The parametric coordinate(s) in which to evaluate
    :param tolerance_in: Knot tolerance for detecting end-point-evaluations
    :param derivatives_in: Number of derivatives to compute

    :param t: evaluation points
    :type  t: [float]
    :return: none
    """
    # wrap everything into c-type datastructures for optimized performence
    cdef unsigned int i,j
    cdef unsigned int n     = len(knots_in)
    cdef double       tol   = tolerance_in
    cdef double[:]    t     = eval_t_in
    cdef double[:]    knots = knots_in
    for j in range(len(t)):
        i = bisect_left(knots, t[j])
        if i < n and abs(knots[i]-t[j]) < tol:
            t[j] = knots[i]
        elif i > 0 and abs(knots[i-1]-t[j]) < tol:
            t[j] = knots[i-1]

