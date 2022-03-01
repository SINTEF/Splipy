from bisect import bisect_right, bisect_left
import numpy as np
cimport numpy as np
import copy
cimport cython


cdef my_bisect_left(np.float_t[:] array, np.float_t value, unsigned int hi):
    cdef unsigned int lo = 0
    cdef unsigned int mid
    while lo < hi:
        mid = (lo + hi) // 2
        if array[mid] < value:
            lo = mid + 1
        else:
            hi = mid
    return lo


cdef my_bisect_right(np.float_t[:] array, np.float_t value, unsigned int hi):
    cdef unsigned int lo = 0
    cdef unsigned int mid
    while lo < hi:
        mid = (lo + hi) // 2
        if value < array[mid]:
            hi = mid
        else:
            lo = mid + 1
    return lo


@cython.boundscheck(False) # turn off bounds-checking for entire function
def evaluate(np.ndarray[np.float_t, ndim=1] knots_in,
             unsigned int p,
             np.ndarray[np.float_t, ndim=1] eval_t_in,
             int periodic,
             np.float_t tol,
             unsigned int d=0,
             bint from_right=True):
    """  Evaluate all basis functions in a given set of points.

    :param knots_in:       Knot vector
    :param p:              Parametrization order (polyonomial degree + 1)
    :param eval_t_in:      The parametric coordinate(s) in which to evaluate
    :param periodic:       Periodicity of basis
    :param tol:            Knot tolerance for detecting end-point-evaluations
    :param d:              Number of derivatives to compute
    :param from_right:     True if evaluation should be done in the limit from above
    :return: Two tuples of all arguments to the scipy sparse csr_matrix
    """

    # wrap everything into c-type datastructures for optimized performance
    cdef np.float_t[:] knots = knots_in
    cdef unsigned int n_all  = len(knots) - p  # number of basis functions (without periodicity)
    cdef unsigned int n      = len(knots) - p - (periodic+1)  # number of basis functions (with periodicity)
    cdef unsigned int m      = len(eval_t_in)
    cdef np.float_t start    = knots[p-1]
    cdef np.float_t end      = knots[n_all]
    cdef np.float_t evalT
    cdef unsigned int mu     = 0

    cdef np.float_t[:] t    = eval_t_in.copy()
    cdef np.float_t[:] data = np.zeros(m*p, dtype=float)
    cdef np.int32_t[:] indices  = np.zeros(m*p, dtype=np.int32)
    cdef np.int32_t[:] indptr   = np.arange(0,m*p+1,p, dtype=np.int32)
    cdef np.float_t[:] M        = np.zeros(p, dtype=float)  # temp storage to keep all the function evaluations
    cdef unsigned int k,q,j,i
    cdef bint right

    if periodic >= 0:
        # Wrap periodic evaluation into domain
        for i in range(len(t)):
            if t[i] < start or t[i] > end:
                t[i] = (t[i] - start) % (end - start) + start
            if abs(t[i] - start) < tol and not from_right:
                t[i] = end
    for i in range(len(t)):
        right = from_right
        evalT = t[i]
        # Special-case the endpoint, so the user doesn't need to
        if abs(t[i] - end) < tol:
            right = False
        # Skip non-periodic evaluation points outside the domain
        if t[i] < start or t[i] > end or (abs(t[i]-start) < tol and not right):
            continue

        # mu = index of last non-zero basis function
        if right:
            mu = my_bisect_right(knots, evalT, n_all+p)
        else:
            mu = my_bisect_left(knots, evalT, n_all+p)
        mu = min(mu, n_all)

        for k in range(p-1):
            M[k] = 0
        M[p-1] = 1  # the last entry is a dummy-zero which is never used
        for q in range(1, p-d):
            j = p-q-1
            k = mu - q -1
            M[j] = M[j] + M[j + 1] * (knots[k + q + 1] - evalT) / (knots[k + q + 1] - knots[k + 1])
            for j in range(p - q , p-1):
                k = mu - p + j  # 'i'-index in global knot vector (ref Hughes book pg.21)
                M[j] = M[j] * (evalT - knots[k]) / (knots[k + q] - knots[k])
                M[j] = M[j] + M[j + 1] * (knots[k + q + 1] - evalT) / (knots[k + q + 1] - knots[k + 1])
            j = p  - 1
            k = mu - 1
            M[j] = M[j] * (evalT - knots[k]) / (knots[k + q] - knots[k])

        for q in range(p-d, p):
            for j in range(p - q - 1, p):
                k = mu - p + j  # 'i'-index in global knot vector (ref Hughes book pg.21)
                if j != p-q-1:
                    M[j] = M[j] * q / (knots[k + q] - knots[k])
                if j != p-1:
                    M[j] = M[j] - M[j + 1] * q / (knots[k + q + 1] - knots[k + 1])


        for j,k in enumerate(range(i*p, (i+1)*p)):
            data[k]    = M[j]
            indices[k] = (mu-p+j) % n
    return (data, indices, indptr), (m,n)


@cython.boundscheck(False) # turn off bounds-checking for entire function
def snap(np.ndarray[np.float_t, ndim=1] knots_in,
         np.ndarray[np.float_t, ndim=1] eval_t_in,
         np.float_t tolerance):
    """  Snap evaluation points to knots if they are sufficiently close
    as given in by state.state.knot_tolerance. This will modify the
    input vector eval_t_in.

    :param knots_in:     Knot vector
    :param eval_t_in:    The parametric coordinate(s) in which to evaluate
    :param tolerance_in: Knot tolerance for detecting end-point-evaluations
    """
    # wrap everything into c-type datastructures for optimized performance
    cdef unsigned int i,j
    cdef unsigned int  n     = len(knots_in)
    cdef np.float_t[:] t     = eval_t_in
    cdef np.float_t[:] knots = knots_in
    for j in range(len(t)):
        i = my_bisect_left(knots, t[j], n)
        if i < n and abs(knots[i]-t[j]) < tolerance:
            t[j] = knots[i]
        elif i > 0 and abs(knots[i-1]-t[j]) < tolerance:
            t[j] = knots[i-1]

