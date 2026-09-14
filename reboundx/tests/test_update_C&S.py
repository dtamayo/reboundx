import numpy as np
import reboundx as rbx
from test_utils import spherical_harmonics_value


def idx(n, m):
    return n * (n + 1) // 2 + m


GM = 6.6743e-11 * 5.9736e24
R_eq = 6378.137e3

N = 4
SIZE = (N + 1) * (N + 2) // 2

# Points spread across n=2..4 (all coefficients defined below are active
# there), chosen away from phi=0, the poles, and lambda multiples of pi/2,
# so cos(m*lambda)/sin(m*lambda) and dP/dphi are all non-trivially exercised.
TEST_POINTS = [
    (np.deg2rad(23.0), 7000e3, np.deg2rad(35.0)),
    (np.deg2rad(-38.0), 8500e3, np.deg2rad(160.0)),
    (np.deg2rad(51.0), 26000e3, np.deg2rad(-70.0)),
]


def _base_coeffs():
    """A handful of active (n,m) terms spanning zonal, sectoral and
    tesseral cases, plus one placeholder coefficient (idx(3,3)) seeded
    with a tiny non-zero value so it stays 'active' and can later be
    turned on for real via update()."""
    C = np.zeros(SIZE)
    S = np.zeros(SIZE)

    C[idx(2, 0)] = -1.0826e-3      # zonal (C_n0, real physical value -- see create() convention)
    C[idx(2, 1)] = 2.5e-7          # tesseral, odd m
    S[idx(2, 1)] = -1.8e-7
    C[idx(2, 2)] = 1.57e-6         # sectoral
    S[idx(2, 2)] = -9.0e-7
    C[idx(4, 0)] = -1.6196e-6      # second zonal, exercises n != 2 too
    C[idx(3, 3)] = 1e-30           # placeholder: "active" but negligible

    return C, S

def model_n_active(C, S):
    """Mirrors create()'s active-coefficient count, for the assertion above."""
    n_active = 0
    n_max = 0
    while (n_max + 1) * (n_max + 2) // 2 < len(C):
        n_max += 1
    for n in range(2, n_max + 1):
        for m in range(0, n + 1):
            k = idx(n, m)
            if k < len(C) and (C[k] != 0.0 or S[k] != 0.0):
                n_active += 1
    return n_active

def test_update_matches_recreate():
    """
    update() must produce EXACTLY the same potential as building a fresh
    model from scratch with the updated coefficients, at every (n,m) that
    was already active at creation time.
    """
    C1, S1 = _base_coeffs()
    model = rbx.spherical_harmonics_model(C1, S1, GM, R_eq)

    # Sanity: capture the "before" potential so we can confirm update()
    # actually changed something (not just silently a no-op).
    V_before = [spherical_harmonics_value(model, phi, r, lam) for phi, r, lam in TEST_POINTS]

    # New coefficient set: every previously-active term gets a new value.
    C2, S2 = _base_coeffs()
    C2[idx(2, 0)] = -1.0827e-3
    C2[idx(2, 1)] = 9.1e-7
    S2[idx(2, 1)] = 3.3e-7
    C2[idx(2, 2)] = -2.2e-6
    S2[idx(2, 2)] = 4.4e-7
    C2[idx(4, 0)] = -1.55e-6
    C2[idx(3, 3)] = 8.2e-7          # now "for real" -- was only a placeholder before

    n_updated = model.update(C2, S2)
    assert n_updated == model_n_active(C1, S1), (
        f"update() reported {n_updated} coefficients updated, expected "
        f"{model_n_active(C1, S1)} (the number active at creation time)."
    )

    model_fresh = rbx.spherical_harmonics_model(C2, S2, GM, R_eq)

    for i, (phi, r, lam) in enumerate(TEST_POINTS):
        V_updated = spherical_harmonics_value(model, phi, r, lam)
        V_fresh = spherical_harmonics_value(model_fresh, phi, r, lam)

        assert V_updated != V_before[i], (
            f"update() had no effect at phi={phi}, r={r}, lam={lam}: "
            f"V is unchanged from before the call."
        )

        rel_err = abs((V_updated - V_fresh) / V_fresh) if V_fresh != 0 else abs(V_updated - V_fresh)
        print(f"phi={np.rad2deg(phi):.1f} r={r:.0f} lam={np.rad2deg(lam):.1f}: "
              f"updated={V_updated:.9e} fresh={V_fresh:.9e} rel_err={rel_err:.2e}")
        assert rel_err < 1e-12, (
            f"update() diverges from a freshly created model at phi={phi}, "
            f"r={r}, lam={lam}: {V_updated} vs {V_fresh}"
        )

    print("update() matches a freshly created model at every test point.")


def test_update_ignores_inactive_coefficients():
    """
    A coefficient that was exactly zero at creation time (so never entered
    model->active) must be silently ignored by update() -- even if it's
    non-zero in the array passed to update(). This documents the
    documented limitation rather than letting it fail silently/confusingly
    if the underlying behaviour ever changes.
    """
    C1, S1 = _base_coeffs()
    assert C1[idx(4, 4)] == 0.0 and S1[idx(4, 4)] == 0.0, "test setup assumption broken"

    model = rbx.spherical_harmonics_model(C1, S1, GM, R_eq)
    phi, r, lam = TEST_POINTS[0]
    V_before = spherical_harmonics_value(model, phi, r, lam)

    C2, S2 = _base_coeffs()
    C2[idx(4, 4)] = 5.0e-7   # was inactive at creation -- update() must ignore this
    S2[idx(4, 4)] = 5.0e-7

    model.update(C2, S2)
    V_after = spherical_harmonics_value(model, phi, r, lam)

    assert V_after == V_before, (
        "update() picked up a coefficient that was inactive at creation "
        "time; it should have been silently ignored. If this is now "
        "intended behaviour, update this test and the docstring in "
        "spherical_harmonics.h together."
    )
    print("update() correctly ignored the newly-nonzero, previously-inactive coefficient.")


def test_update_zeroing_a_coefficient_is_not_special_cased():
    """
    A coefficient that WAS active at creation can be updated down to
    exactly 0.0 without any special handling -- it should just contribute
    zero force/potential from then on, matching a fresh model built with
    that coefficient at 0.0.
    """
    C1, S1 = _base_coeffs()
    model = rbx.spherical_harmonics_model(C1, S1, GM, R_eq)

    C2, S2 = _base_coeffs()
    C2[idx(2, 1)] = 0.0
    S2[idx(2, 1)] = 0.0

    model.update(C2, S2)
    model_fresh = rbx.spherical_harmonics_model(C2, S2, GM, R_eq)

    for phi, r, lam in TEST_POINTS:
        V_updated = spherical_harmonics_value(model, phi, r, lam)
        V_fresh = spherical_harmonics_value(model_fresh, phi, r, lam)
        rel_err = abs((V_updated - V_fresh) / V_fresh) if V_fresh != 0 else abs(V_updated - V_fresh)
        assert rel_err < 1e-12, (
            f"Zeroing an active coefficient via update() diverges from a "
            f"fresh model at phi={phi}, r={r}, lam={lam}: {V_updated} vs {V_fresh}"
        )
    print("Zeroing an active coefficient via update() matches a fresh model.")




if __name__ == "__main__":
    test_update_matches_recreate()
    test_update_ignores_inactive_coefficients()
    test_update_zeroing_a_coefficient_is_not_special_cased()
