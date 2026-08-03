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

        if self._freed or self._ptr is None:
            raise RuntimeError("Can not attach a spherical_harmonics_model already freed")

        particle.params["spherical_harmonics_model"] = self._ptr

        sim = particle.sim
        if sim is None:
            raise ValueError("The particle has to be part of a simulation before the attach()")

        if not hasattr(sim, "_rebx_anchors"):
            sim._rebx_anchors = []
        sim._rebx_anchors.append(self)

        self._attached_to = particle