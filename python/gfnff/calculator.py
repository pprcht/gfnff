# This file is part of gfnff.
# SPDX-Identifier: LGPL-3.0-or-later
"""
Low-level Python wrapper around the GFN-FF C API.

All quantities use the native C API units:
  - positions / lattice : Bohr
  - energy              : Hartree
  - gradient            : Eh / Bohr
"""

import ctypes

import numpy as np

from ._lib import _CGFNFFCalculator, _lib


class GFNFFCalculator:
    """Wraps the GFN-FF C API for a single molecular or periodic system.

    Parameters
    ----------
    numbers:
        Atomic numbers, shape ``(nat,)``.
    positions_bohr:
        Cartesian coordinates in Bohr, shape ``(nat, 3)``, C-contiguous.
    charge:
        Total molecular charge (integer).
    printlevel:
        Verbosity of Fortran output (0 = silent, 1 = errors only,
        2 = informational, 3 = verbose).
    solvent:
        Implicit solvent name (e.g. ``"h2o"``, ``"acetone"``).
        Empty string or ``None`` disables solvation.
    lattice:
        Lattice vectors in Bohr, shape ``(3, 3)``, C-contiguous.
        Each *row* is one lattice vector (same convention as the C API).
        ``None`` for non-periodic systems.
    npbc:
        Number of periodic dimensions (0–3).  Ignored when
        ``lattice`` is ``None``.
    """

    def __init__(
        self,
        numbers,
        positions_bohr,
        *,
        charge: int = 0,
        printlevel: int = 0,
        solvent: str = "",
        lattice=None,
        npbc: int = 0,
    ):
        at = np.asarray(numbers, dtype=np.int32)
        xyz = np.ascontiguousarray(positions_bohr, dtype=np.float64)
        if xyz.ndim != 2 or xyz.shape[1] != 3:
            raise ValueError(
                f"positions_bohr must have shape (nat, 3), got {xyz.shape}"
            )
        if at.shape[0] != xyz.shape[0]:
            raise ValueError("numbers and positions_bohr must have the same length")

        nat = ctypes.c_int(len(at))
        solvent_bytes = (solvent or "").encode()

        if lattice is not None:
            lat = np.ascontiguousarray(lattice, dtype=np.float64)
            if lat.shape != (3, 3):
                raise ValueError(f"lattice must have shape (3, 3), got {lat.shape}")
            self._handle = _lib.c_gfnff_calculator_init_pbc(
                nat, at, xyz,
                ctypes.c_int(charge),
                ctypes.c_int(printlevel),
                lat,
                ctypes.c_int(npbc),
            )
        else:
            self._handle = _lib.c_gfnff_calculator_init(
                nat, at, xyz,
                ctypes.c_int(charge),
                ctypes.c_int(printlevel),
                solvent_bytes,
            )

        if self._handle.ptr is None:
            raise RuntimeError(
                "GFN-FF initialization failed (c_gfnff_calculator_init returned NULL). "
                "Check printlevel > 0 for error details."
            )
        self._alive = True

    # ------------------------------------------------------------------

    def singlepoint(self, numbers, positions_bohr, lattice=None):
        """Compute energy, gradient and stress tensor for the given geometry.

        Parameters
        ----------
        numbers:
            Atomic numbers, shape ``(nat,)``.
        positions_bohr:
            Coordinates in Bohr, shape ``(nat, 3)``.
        lattice:
            Lattice vectors in Bohr, shape ``(3, 3)``.  Must be provided
            for periodic systems (same ``npbc`` as at initialization).

        Returns
        -------
        energy : float
            Total energy in Hartree.
        gradient : ndarray, shape ``(nat, 3)``
            Energy gradient in Eh / Bohr.
        sigma : ndarray, shape ``(3, 3)``
            Stress tensor in Hartree.  Zero for non-periodic systems.
        """
        self._check_alive()

        at = np.asarray(numbers, dtype=np.int32)
        xyz = np.ascontiguousarray(positions_bohr, dtype=np.float64)
        nat = len(at)

        energy_out = ctypes.c_double(0.0)
        gradient = np.zeros((nat, 3), dtype=np.float64, order="C")
        sigma = np.zeros((3, 3), dtype=np.float64, order="C")
        iostat = ctypes.c_int(0)

        _lib.c_gfnff_calculator_singlepoint(
            ctypes.byref(self._handle),
            ctypes.c_int(nat),
            at,
            xyz,
            ctypes.byref(energy_out),
            gradient,
            sigma,
            ctypes.byref(iostat),
        )

        if iostat.value != 0:
            raise RuntimeError(
                f"GFN-FF singlepoint failed with iostat = {iostat.value}"
            )

        return float(energy_out.value), gradient, sigma

    # ------------------------------------------------------------------

    def print_results(self, iunit: int = 6):
        """Print the GFN-FF energy decomposition to a Fortran unit.

        Parameters
        ----------
        iunit:
            Fortran I/O unit number.  Use ``6`` for standard output.
        """
        self._check_alive()
        _lib.c_gfnff_calculator_results(
            ctypes.byref(self._handle),
            ctypes.c_int(iunit),
        )

    # ------------------------------------------------------------------

    def _check_alive(self):
        if not self._alive or self._handle.ptr is None:
            raise RuntimeError("GFNFFCalculator has already been deallocated.")

    def deallocate(self):
        """Explicitly free Fortran-side memory."""
        if getattr(self, "_alive", False) and self._handle.ptr is not None:
            _lib.c_gfnff_calculator_deallocate(ctypes.byref(self._handle))
            self._alive = False

    def __del__(self):
        self.deallocate()

    def __enter__(self):
        return self

    def __exit__(self, *_):
        self.deallocate()
