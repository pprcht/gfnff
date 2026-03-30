# This file is part of gfnff.
# SPDX-Identifier: LGPL-3.0-or-later
"""
ctypes loader and raw C API prototypes for libgfnff.

All five public functions from gfnff_interface_c.h are wrapped here.
Units match the C API exactly: Bohr, Hartree, Eh/Bohr.
"""

import ctypes
import ctypes.util
import pathlib

import numpy as np
import numpy.ctypeslib as npct


# ---------------------------------------------------------------------------
# Library loading
# ---------------------------------------------------------------------------

def _load_library() -> ctypes.CDLL:
    # 1. Wheel / editable install: .so lives next to this file
    pkg_dir = pathlib.Path(__file__).parent
    candidates = (
        list(pkg_dir.glob("libgfnff.so*"))   # Linux: libgfnff.so.0.0.1
        + list(pkg_dir.glob("libgfnff*.dylib"))  # macOS: libgfnff.0.0.1.dylib
        + list(pkg_dir.glob("gfnff.dll"))
    )
    if candidates:
        return ctypes.CDLL(str(candidates[0]))

    # 2. LD_LIBRARY_PATH / system library
    name = ctypes.util.find_library("gfnff")
    if name:
        return ctypes.CDLL(name)

    raise OSError(
        "libgfnff shared library not found. "
        "Install the gfnff package or set LD_LIBRARY_PATH to the directory "
        "containing libgfnff.so."
    )


_lib = _load_library()


# ---------------------------------------------------------------------------
# Opaque calculator struct
# ---------------------------------------------------------------------------

class _CGFNFFCalculator(ctypes.Structure):
    """Mirror of ``typedef struct { void *ptr; } c_gfnff_calculator;``"""
    _fields_ = [("ptr", ctypes.c_void_p)]


# ---------------------------------------------------------------------------
# Convenience ndpointer types
# ---------------------------------------------------------------------------

# double xyz[nat][3]  or  double gradient[nat][3]  -- shape checked at runtime
_xyz_type = npct.ndpointer(dtype=np.float64, ndim=2, flags="C_CONTIGUOUS")

# int at[nat]
_at_type = npct.ndpointer(dtype=np.int32, ndim=1, flags="C_CONTIGUOUS")

# double lattice[3][3]
_lattice_type = npct.ndpointer(dtype=np.float64, ndim=2, flags="C_CONTIGUOUS")


# ---------------------------------------------------------------------------
# Function prototypes
# ---------------------------------------------------------------------------

# c_gfnff_calculator c_gfnff_calculator_init(
#     int nat, int *at, double (*xyz)[3],
#     int ichrg, int printlevel, const char *solvent);
_lib.c_gfnff_calculator_init.restype = _CGFNFFCalculator
_lib.c_gfnff_calculator_init.argtypes = [
    ctypes.c_int,    # nat
    _at_type,        # at[nat]
    _xyz_type,       # xyz[nat][3]
    ctypes.c_int,    # ichrg
    ctypes.c_int,    # printlevel
    ctypes.c_char_p, # solvent (null-terminated)
]

# c_gfnff_calculator c_gfnff_calculator_init_pbc(
#     int nat, int *at, double (*xyz)[3],
#     int ichrg, int printlevel,
#     double lattice[3][3], int npbc);
_lib.c_gfnff_calculator_init_pbc.restype = _CGFNFFCalculator
_lib.c_gfnff_calculator_init_pbc.argtypes = [
    ctypes.c_int,    # nat
    _at_type,        # at[nat]
    _xyz_type,       # xyz[nat][3]
    ctypes.c_int,    # ichrg
    ctypes.c_int,    # printlevel
    _lattice_type,   # lattice[3][3]
    ctypes.c_int,    # npbc
]

# void c_gfnff_calculator_singlepoint(
#     c_gfnff_calculator *calculator,
#     int nat, int *at, double (*xyz)[3],
#     double *energy, double (*gradient)[3],
#     double sigma[3][3], int *iostat);
_lib.c_gfnff_calculator_singlepoint.restype = None
_lib.c_gfnff_calculator_singlepoint.argtypes = [
    ctypes.POINTER(_CGFNFFCalculator),  # calculator
    ctypes.c_int,                        # nat
    _at_type,                            # at[nat]
    _xyz_type,                           # xyz[nat][3]
    ctypes.POINTER(ctypes.c_double),     # energy (out)
    _xyz_type,                           # gradient[nat][3] (out)
    _lattice_type,                       # sigma[3][3] (out); zero for non-PBC
    ctypes.POINTER(ctypes.c_int),        # iostat (out)
]

# void c_gfnff_calculator_results(c_gfnff_calculator *calculator, int iunit);
_lib.c_gfnff_calculator_results.restype = None
_lib.c_gfnff_calculator_results.argtypes = [
    ctypes.POINTER(_CGFNFFCalculator),  # calculator
    ctypes.c_int,                        # iunit (6 = Fortran stdout)
]

# void c_gfnff_calculator_deallocate(c_gfnff_calculator *calculator);
_lib.c_gfnff_calculator_deallocate.restype = None
_lib.c_gfnff_calculator_deallocate.argtypes = [
    ctypes.POINTER(_CGFNFFCalculator),  # calculator
]
