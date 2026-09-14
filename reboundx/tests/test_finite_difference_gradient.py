import numpy as np
import reboundx as rbx
from reboundx.tests.test_potential import independent_potential
from test_utils import spherical_harmonics_value, spherical_harmonics_acceleration
from math import factorial

def idx(n, m):
    return n * (n + 1) // 2 + m


def test_c21_matches_independent_reference():
    """
    This test compares the model's potential directly against an independent 
    scipy-based reference, for an odd m coefficient (m=1), where (-1)^m = -1 
    makes a phase convention mismatch visible.
    """
    N = 4
    size = (N + 1) * (N + 2) // 2
    GM = 6.6743e-11 * 5.9736e24
    R_eq = 6378.137e3

    C = np.zeros(size)
    S = np.zeros(size)
    C[idx(2, 1)] = 2.5e-7
    S[idx(2, 1)] = -1.8e-7

    model = rbx.spherical_harmonics_model(C, S, GM, R_eq)

    test_points = [
        (np.deg2rad(23.0), 7000e3, np.deg2rad(35.0)),
        (np.deg2rad(-38.0), 8500e3, np.deg2rad(160.0)),
        (np.deg2rad(51.0), 26000e3, np.deg2rad(-70.0)),
    ]

    for phi, r, lam in test_points:
        V_model = spherical_harmonics_value(model, phi, r, lam)
        V_ref = independent_potential(r, phi, lam, N, C, S, GM, R_eq)
        rel_err = abs((V_model - V_ref) / V_ref)
        print(f"phi={np.rad2deg(phi):.1f} r={r:.0f}: model={V_model:.6e} ref={V_ref:.6e} rel_err={rel_err:.2e}")
        assert rel_err < 1e-9, f"CS phase mismatch at phi={phi}, r={r}: {V_model} vs {V_ref}"

    print("C21/S21 potential matches independent (scipy, CS-phase-included) reference.")

def test_c21_is_normalized_matches_unnormalized():
    """
    Same physical C21/S21 as test_c21_matches_independent_reference, but
    passed through is_normalized=True instead. The two must give the
    IDENTICAL potential, since they're mathematically the same coefficient
    expressed in two different conventions: Cbar_nm = C_nm / N_nm(n,m).
    """
    N = 4
    size = (N + 1) * (N + 2) // 2
    GM = 6.6743e-11 * 5.9736e24
    R_eq = 6378.137e3

    C21_raw = 2.5e-7
    S21_raw = -1.8e-7

    norm_factor_21 = N_nm(2, 1)   
    C21_bar = C21_raw / norm_factor_21
    S21_bar = S21_raw / norm_factor_21

    C_unnorm = np.zeros(size); S_unnorm = np.zeros(size)
    C_unnorm[idx(2, 1)] = C21_raw
    S_unnorm[idx(2, 1)] = S21_raw

    C_norm = np.zeros(size); S_norm = np.zeros(size)
    C_norm[idx(2, 1)] = C21_bar
    S_norm[idx(2, 1)] = S21_bar

    model_unnorm = rbx.spherical_harmonics_model(C_unnorm, S_unnorm, GM, R_eq, is_normalized=False)
    model_norm = rbx.spherical_harmonics_model(C_norm, S_norm, GM, R_eq, is_normalized=True)

    test_points = [
        (np.deg2rad(23.0), 7000e3, np.deg2rad(35.0)),
        (np.deg2rad(-38.0), 8500e3, np.deg2rad(160.0)),
        (np.deg2rad(51.0), 26000e3, np.deg2rad(-70.0)),
    ]

    for phi, r, lam in test_points:
        V_unnorm = spherical_harmonics_value(model_unnorm, phi, r, lam)
        V_norm = spherical_harmonics_value(model_norm, phi, r, lam)
        rel_err = abs((V_norm - V_unnorm) / V_unnorm)
        print(f"phi={np.rad2deg(phi):.1f} r={r:.0f}: unnorm={V_unnorm:.6e} norm={V_norm:.6e} rel_err={rel_err:.2e}")
        assert rel_err < 1e-12, (
            f"is_normalized=True diverges from is_normalized=False (should be "
            f"the exact same physics, two conventions) at phi={phi}, r={r}: "
            f"{V_norm} vs {V_unnorm}"
        )

    print("is_normalized=True matches is_normalized=False for the equivalent C21/S21.")


def N_nm(n, m):
    """ Normalization factor, same formula as N_nm() in spherical_harmonics.c """
    delta_0m = 1 if m == 0 else 0
    return np.sqrt((2 - delta_0m) * (2 * n + 1) * factorial(n - m) / factorial(n + m))


if __name__ == "__main__":
    test_c21_matches_independent_reference()
    test_c21_is_normalized_matches_unnormalized()