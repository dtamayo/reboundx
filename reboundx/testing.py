"""
These are structures anad functions mainly used for testing from the Python side and for use in unit tests.
"""

from . import clibreboundx
from ctypes import c_void_p, c_int, c_long, Structure, byref, c_char_p
import warnings
from reboundx.extras import REBX_BINARY_WARNINGS

REBX_BINARY_FIELD_TYPE = {
        0: 'None',
        1: 'Operator',
        2: 'Particle',
        3: 'Rebx Structure',
        4: 'Param',
        5: 'Name',
        6: 'Param type',
        7: 'Param valaue',
        8: 'END',
        9: 'Particle index',
        10: 'Rebx integraator',
        11: 'Force type',
        12: 'Operator type',
        13: 'Step',
        14: 'Step dt fraction',
        15: 'Registered Param',
        16: 'Additional force',
        17: 'Param list',
        18: 'Registered Params',
        19: 'Allocated forces',
        20: 'Allocated operators',
        21: 'Additional forces',
        22: 'Pre timestep modifications',
        23: 'Post timestep modifications',
        24: 'Particles',
        25: 'Force',
        26: 'Snapshot',
        }

class BinaryField(Structure):
    _fields_ =  [  ("_type", c_int),
                    ("size", c_long)]
    @property
    def type(self):
        return REBX_BINARY_FIELD_TYPE[self._type]

    def __repr__(self):
        return 'Type: {0}, Size: {1}'.format(self.type, self.size)

def inspect_binary(filename):
    w = c_int(0)
    clibreboundx.rebx_input_inspect_binary.restype = c_void_p
    inf = clibreboundx.rebx_input_inspect_binary(c_char_p(filename.encode('ascii')), byref(w))
    for majorerror, value, message in REBX_BINARY_WARNINGS:
        if w.value & value:
            if majorerror:
                raise RuntimeError(message)
            else:
                warnings.warn(message, RuntimeWarning)
    return inf

def read_binary_field(inf):
    clibreboundx.rebx_input_read_binary_field.restype = BinaryField
    return clibreboundx.rebx_input_read_binary_field(c_void_p(inf))

def skip_binary_field(inf, size):
    clibreboundx.rebx_input_skip_binary_field(c_void_p(inf), size)

