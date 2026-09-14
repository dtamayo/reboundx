from ctypes import c_int, c_double
import numpy as np
from . import clibreboundx


class spherical_harmonics_model(object):
    """
    Spherical Harmonics Model for REBOUNDx.

    The user only interacts with this class - internally returns 
    the calls of rebx_spherical_harmonics_create / rebx_spherical_harmonics_free /
    rebx_set_param_pointer via ctypes.

    Parameters
    ----------
    C, S : list or float array
    GM : float
    R_eq : float
    is_normalized : bool, optional

    Example
    -------
    >>> model = rbx.spherical_harmonics_model(C, S, GM, R_eq)
    >>> model.attach(central_particle)
    """

    def __init__(self, C, S, GM, R_eq, is_normalized=False):
        
        self._ptr = None
        self._freed = True

        size = len(C)
        N = int(round((-3 + np.sqrt(1 + 8 * size)) / 2))

        if (N + 1) * (N + 2) // 2 != size:
            raise ValueError(
                f"Invalid coefficient array length ({size}). "
                "Expected (N+1)(N+2)/2 elements."
            )
        
        if len(C) != len(S):
            raise ValueError("C and S must have the same length.")

        C_arr = (c_double * size)()
        S_arr = (c_double * size)()

        for i in range(len(C)):
            C_arr[i] = C[i]
            S_arr[i] = S[i]

        self._ptr = clibreboundx.rebx_spherical_harmonics_create(
            c_int(N), C_arr, S_arr, c_double(GM), c_double(R_eq), c_int(1 if is_normalized else 0)
        )

        if not self._ptr:
            raise RuntimeError("rebx_spherical_harmonics_create returned NULL")

        self._freed = False
        self._size = size   # NEW: To validate update()
        self._N = N 

    def update(self, C, S):
        """
        Re-read C and S at every (n,m) coefficient that was already active
        at creation time, and update the model's internal coefficients
        in place -- without recreating the model.

        C and S must have the SAME length and dense idx(n,m) layout, and
        the SAME convention (raw C_nm/S_nm, is_normalized) as originally
        passed to the constructor. Coefficients that were zero at creation
        time are silently ignored, even if non-zero here -- see the C
        header docstring for rebx_spherical_harmonics_update().

        Returns
        -------
        int
            Number of coefficients updated.
        """
        if self._freed or self._ptr is None:
            raise RuntimeError("Cannot update a spherical_harmonics_model that has been freed")

        if len(C) != self._size or len(S) != self._size:
            raise ValueError(
                f"C and S must have length {self._size} (same as at creation), "
                f"got {len(C)} and {len(S)}."
            )

        C_arr = (c_double * self._size)()
        S_arr = (c_double * self._size)()
        for i in range(self._size):
            C_arr[i] = C[i]
            S_arr[i] = S[i]

        n_updated = clibreboundx.rebx_spherical_harmonics_update(self._ptr, C_arr, S_arr)
        return n_updated


    def __del__(self):
        if getattr(self, "_freed", True):
            return
        if getattr(self, "_ptr", None):
            try:
                clibreboundx.rebx_spherical_harmonics_free(self._ptr)
            except Exception:
                pass
            self._freed = True

def attach(self, particle):
    """
    Attach this model to a particle by writing its C pointer into
    particle.params["spherical_harmonics_model"].

    LIFETIME WARNING: this method does NOT protect the model from Python's
    garbage collector. It only writes the pointer into particle.params; it
    keeps no Python reference to `self` anywhere. If nothing else in your
    script holds a reference to this model object, Python is free to
    garbage-collect it (and free the underlying C memory) even while the
    simulation keeps running with a now-dangling pointer.

    You are responsible for keeping this model -- and the force object
    returned by rebx.load_force(...) -- alive for as long as the
    simulation runs. The established REBOUNDx pattern for this is to hang
    them off the Extras instance under a name starting with "_" (the
    convention that keeps such references exempt from REBOUND's attribute
    validation):

        rebx = reboundx.Extras(sim)
        gh = rebx.load_force("gravitational_harmonics")
        model = spherical_harmonics_model(C, S, GM, R_eq)
        model.attach(sim.particles["Earth"])
        rebx._gh_force = gh                        # keeps the force alive
        rebx._spherical_harmonics_model = model    # keeps the model alive
        rebx.add_force(gh)

    Example
    -------
    >>> model = rbx.spherical_harmonics_model(C, S, GM, R_eq)
    >>> model.attach(sim.particles["Earth"])
    >>> rebx._spherical_harmonics_model = model  # <- required, see LIFETIME WARNING above
    """
    if self._freed or self._ptr is None:
        raise RuntimeError("Cannot attach a spherical_harmonics_model that has already been freed")

    if particle.sim is None:
        raise ValueError("The particle must be part of a simulation before calling attach()")

    particle.params["spherical_harmonics_model"] = self._ptr
    self._attached_to = particle