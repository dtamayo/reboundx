import numpy as np
from scipy.special import assoc_legendre_p_all
from ctypes import c_int, c_double, POINTER, byref
from reboundx import clibreboundx


def idx(n, m):
    return n * (n + 1) // 2 + m


def factorial(n):
    f = 1.0
    for i in range(1, n + 1):
        f *= i
    return f


def N_nm(n, m):
    delta = 1 if m == 0 else 0
    return np.sqrt((2 - delta) * (2 * n + 1) * factorial(n - m) / factorial(n + m))


def reference_P(N, phi):
    """ Fully normalized, Condon-Shortley-phase-included associated Legendre
    functions, matching the ICGEM/geodesy convention that EGM2008 and
    friends are published in. Independent of anything in legendre.c. """
    x = np.sin(phi)
    P_all = assoc_legendre_p_all(N, N, x, norm=False)[0]
    P = np.zeros(idx(N, N) + 1)
    for n in range(0, N + 1):
        for m in range(0, n + 1):
            P[idx(n, m)] = N_nm(n, m) * P_all[n, m]
    return P


# --- ctypes bindings straight to the C recursion, no potential/acceleration
# wrapper in between, so a bug here can't hide behind self-consistency with
# any other part of the pipeline. ---

clibreboundx.compute_recursion_coef.argtypes = [
    c_int, POINTER(c_double), POINTER(c_double),
    POINTER(c_double), POINTER(c_double), POINTER(c_double),
]
clibreboundx.compute_recursion_coef.restype = None

clibreboundx.compute_normalized_legendre.argtypes = [
    c_double, c_int, POINTER(c_double), POINTER(c_double),
    POINTER(c_double), POINTER(c_double), POINTER(c_double),
]
clibreboundx.compute_normalized_legendre.restype = None

clibreboundx.compute_normalized_legendre_derivative.argtypes = [
    c_double, c_int, POINTER(c_double), POINTER(c_double),
    POINTER(c_double), POINTER(c_double), POINTER(c_double),
    POINTER(c_double), POINTER(c_double),
]
clibreboundx.compute_normalized_legendre_derivative.restype = None


def _make_coef_arrays(N):
    size = idx(N, N) + 1
    A = (c_double * size)()
    B = (c_double * size)()
    D = (c_double * (N + 1))()
    E = (c_double * (N + 1))()
    F = (c_double * size)()
    clibreboundx.compute_recursion_coef(c_int(N), A, B, D, E, F)
    return A, B, D, E, F


def _compute_P(xp, N, A, B, D, E):
    size = idx(N, N) + 1
    P = (c_double * size)()
    clibreboundx.compute_normalized_legendre(c_double(xp), c_int(N), A, B, D, E, P)
    return np.array(P[:])


def _compute_dP(phi, N, A, B, D, E, F, P_ctypes):
    size = idx(N, N) + 1
    dP = (c_double * size)()
    clibreboundx.compute_normalized_legendre_derivative(
        c_double(phi), c_int(N), A, B, D, E, F, P_ctypes, dP
    )
    return np.array(dP[:])


N = 8
TEST_LATITUDES_DEG = [23.0, -38.0, 51.0, 5.0, -89.0, 60.0]
TOL_P = 1e-10
TOL_DP = 1e-6   # finite-difference floor, not a formula-precision floor -- see note below


def test_legendre_P_matches_scipy_reference():
    """
    compute_normalized_legendre() must match a fully normalized, Condon-
    Shortley-phase-included associated Legendre reference (scipy) across
    the full (n,m) triangle, not just m==0 (which test_potential.py already
    covers) or a single isolated (n,m) (which the C21/S21 tests cover).

    This is the test that would have caught the D[1] bug (missing sqrt(2)
    correction for the m=0 -> m=1 delta_0m normalization jump) and the
    missing Condon-Shortley sign on the diagonal recursion.
    """
    A, B, D, E, F = _make_coef_arrays(N)

    for phi_deg in TEST_LATITUDES_DEG:
        phi = np.deg2rad(phi_deg)
        P_code = _compute_P(np.sin(phi), N, A, B, D, E)
        P_ref = reference_P(N, phi)

        for n in range(N + 1):
            for m in range(n + 1):
                k = idx(n, m)
                if abs(P_ref[k]) < 1e-12:
                    continue
                rel_err = abs((P_code[k] - P_ref[k]) / P_ref[k])
                assert rel_err < TOL_P, (
                    f"P mismatch at phi={phi_deg} deg, (n,m)=({n},{m}): "
                    f"code={P_code[k]} vs scipy={P_ref[k]} (rel_err={rel_err:.2e})"
                )
    print(f"compute_normalized_legendre() matches scipy for all (n,m) up to N={N}, "
          f"at {len(TEST_LATITUDES_DEG)} latitudes.")


def test_legendre_derivative_matches_finite_difference_of_reference():
    """
    compute_normalized_legendre_derivative() must match the finite-
    difference derivative of the EXTERNAL reference P (not of the model's
    own P) -- this is what makes the test able to catch a bug that
    affects P and dP consistently, since comparing dP against a
    finite-difference of the model's OWN P can pass even if both are
    wrong in the same way.

    Points very close to phi=0 are excluded: finite-difference truncation
    noise there dominates any real formula error (verified independently
    with mpmath arbitrary-precision arithmetic; see the conversation this
    test originated from). This is a numerical floor of the validation
    method, not a floor on the formula's actual accuracy.
    """
    A, B, D, E, F = _make_coef_arrays(N)
    h = 1e-6

    for phi_deg in TEST_LATITUDES_DEG:
        phi = np.deg2rad(phi_deg)

        P_ref = reference_P(N, phi)
        P_ref_minus = reference_P(N, phi - h)
        P_ref_plus = reference_P(N, phi + h)
        dP_fd = (P_ref_plus - P_ref_minus) / (2 * h)

        P_ref_ctypes = (c_double * len(P_ref))(*P_ref)
        dP_code = _compute_dP(phi, N, A, B, D, E, F, P_ref_ctypes)

        for n in range(N + 1):
            for m in range(n + 1):
                k = idx(n, m)
                if abs(dP_fd[k]) < 1e-6:
                    continue  # too close to a zero-crossing for a meaningful relative check
                rel_err = abs((dP_code[k] - dP_fd[k]) / dP_fd[k])
                assert rel_err < TOL_DP, (
                    f"dP mismatch at phi={phi_deg} deg, (n,m)=({n},{m}): "
                    f"code={dP_code[k]} vs fd(scipy)={dP_fd[k]} (rel_err={rel_err:.2e})"
                )
    print(f"compute_normalized_legendre_derivative() matches finite-difference "
          f"of scipy's reference for all (n,m) up to N={N}, "
          f"at {len(TEST_LATITUDES_DEG)} latitudes.")


if __name__ == "__main__":
    test_legendre_P_matches_scipy_reference()
    test_legendre_derivative_matches_finite_difference_of_reference()