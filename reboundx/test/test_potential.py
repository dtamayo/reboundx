import numpy as np
from scipy.special import assoc_legendre_p_all, factorial
import reboundx as rbx
from test_utils import spherical_harmonics_value


def idx(n, m):
    return n * (n + 1) // 2 + m


def independent_potential(r, phi, lam, N, C, S, GM, R_eq):

    x = np.sin(phi)
    V = 0.0
    P_all = assoc_legendre_p_all(N, N, x, norm=False)

    assert P_all.shape == (1, N + 1, 2 * N + 1), f"Non expected shape: {P_all.shape}"

    for n in range(2, N + 1):
        for m in range(0, n + 1):
            Pnm_raw = P_all[0, n, m]
            k = 2 if m > 0 else 1
            norm = np.sqrt((2 * n + 1) * k * factorial(n - m) / factorial(n + m))
            Pnm = norm * Pnm_raw
            i = idx(n, m)
            V += (R_eq / r) ** n * Pnm * (C[i] * np.cos(m * lam) + S[i] * np.sin(m * lam))
    return GM / r * V


def test_potential_matches_independent_reference():
    N = 15
    size = (N + 1) * (N + 2) // 2
    GM = 6.6743e-11 * 5.9736e24
    R_eq = 6378.137e3

    J2 = 1.0826e-3
    J4 = -1.6196e-6

    C = np.zeros(size, dtype=np.float64)
    S = np.zeros(size, dtype=np.float64)
    C[idx(2, 0)] = -J2 / np.sqrt(5)
    C[idx(4, 0)] = -J4 / 3

    model = rbx.spherical_harmonics_model(C, S, GM, R_eq)

    test_points = [
        (np.deg2rad(20), 7000e3, np.deg2rad(10)),
        (np.deg2rad(-45), 8500e3, np.deg2rad(200)),
        (np.deg2rad(60), 42164e3, np.deg2rad(-30)),  
    ]

    for phi, r, lam in test_points:
        V_model = spherical_harmonics_value(model, phi, r, lam)
        V_ref = independent_potential(r, phi, lam, N, C, S, GM, R_eq)
        rel_err = abs((V_model - V_ref) / V_ref)
        print(f"phi={np.rad2deg(phi):.1f} r={r:.0f}: model={V_model:.6e} ref={V_ref:.6e} rel_err={rel_err:.2e}")
        assert rel_err < 1e-9, f"Divergence in phi={phi}, r={r}: {V_model} vs {V_ref}"

    print("Model's Potencial is the same as the independent reference.")


if __name__ == "__main__":
    test_potential_matches_independent_reference()