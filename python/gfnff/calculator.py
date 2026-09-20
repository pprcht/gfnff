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
import os
from enum import IntEnum

import numpy as np

from ._lib import (
    _CGFNFFCalculator,
    _lib,
    NO_CHARGES,
    PBC_UNSUPPORTED,
    VERSION_DEFAULT,
    toml_available,
)


#: Names of the parametrisation versions the library knows about.
_VERSION_NAMES = (
    "angewChem2020",
    "angewChem2020_1",
    "angewChem2020_2",
    "mcgfnff2023",
    "conformer2020",
    "harmonic2020",
)

#: Values are read back from the shared library rather than written out here,
#: so this enum cannot drift from TGFFVersionEnum in the Fortran source.
_version_values = {
    name: _lib.c_gfnff_version_from_name(name.encode()) for name in _VERSION_NAMES
}
_unknown = [name for name, val in _version_values.items() if val == VERSION_DEFAULT]
if _unknown:
    raise ImportError(
        f"libgfnff does not recognise the parametrisation version(s) {_unknown}. "
        "The Python bindings and the shared library are out of step."
    )

Version = IntEnum("Version", {"default": VERSION_DEFAULT, **_version_values})
Version.__doc__ = """GFN-FF parametrisation versions.

``default`` leaves the choice to the library (currently ``angewChem2020_2``).
``mcgfnff2023`` is for molecular crystals and is rejected for non-periodic
systems.  ``harmonic2020`` replaces the bond potential with a harmonic one.
``conformer2020`` reproduces ``angewChem2020_2`` wherever the bonds are near
their equilibrium lengths, but its bonds cannot dissociate.
"""


def _resolve_version(version) -> int:
    """Accept a ``Version``, a plain int, a name, or ``None``."""
    if version is None:
        return VERSION_DEFAULT
    if isinstance(version, str):
        try:
            return int(Version[version])
        except KeyError:
            raise ValueError(
                f"unknown GFN-FF version {version!r}; "
                f"expected one of {', '.join(_VERSION_NAMES)}"
            ) from None
    return int(version)


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
    accuracy:
        Cutoff/precision factor.  Larger is looser and faster; above ``1.0``
        the EEQ system is solved in single precision.  ``None`` uses the
        library default (``0.1``, or ``2.0`` above 10000 atoms).  Must be
        given here rather than changed later, because it also sets the
        thresholds used to build the topology.
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
        fragments=None,
        ref_charges=None,
        accuracy: float | None = None,
        version=None,
        parametrisation: str | None = None,
        bond_matrix=None,
    ):
        """
        version:
            Parametrisation version: a ``Version`` member, its integer value,
            or its name (e.g. ``"mcgfnff2023"``).  ``None`` uses the library
            default.  ``mcgfnff2023`` is rejected for non-periodic systems.
        parametrisation:
            Path to a parameter file.  ``None`` uses the internal set for
            ``version``.  A ``.toml`` path is applied as an *overlay* on that
            set, so the file need only name the keys it changes; any other
            path is read as the legacy flat format.  TOML support requires a
            build with toml-f.
        fragments:
            Optional per-atom fragment index, shape ``(nat,)`` (int). When given,
            GFN-FF forms no bonds between atoms of differing fragments, enforcing
            a host-defined fragmentation instead of the automatic detection.
        ref_charges:
            Optional per-atom reference charges, shape ``(nat,)`` (float). Summed
            over each fragment to set the per-fragment EEQ charge constraint.
            No charge model is invoked; the caller supplies the array directly.
        bond_matrix:
            Optional molecular graph, shape ``(nat, nat)`` of integer bond
            orders.  A nonzero element declares a bond.  When given it replaces
            GFN-FF's distance-based bond perception entirely, so the
            connectivity need not be supported by the coordinates -- which is
            the point: it lets a structure be built from a graph plus an
            arbitrary starting geometry.  Must be symmetric with a zero
            diagonal, non-negative, at most 41 bonds per atom, and
            non-periodic.  Only the zero/nonzero pattern is used today; the
            orders are validated and stored but have no consumer, since GFN-FF
            derives its own pi bond orders.
        """
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

        # validate and keep references to optional bundle arrays alive
        frag_arr = None
        frag_ptr = None
        if fragments is not None:
            frag_arr = np.ascontiguousarray(fragments, dtype=np.int32)
            if frag_arr.shape != (len(at),):
                raise ValueError(
                    f"fragments must have shape ({len(at)},), got {frag_arr.shape}"
                )
            frag_ptr = frag_arr.ctypes.data_as(ctypes.c_void_p)

        refq_arr = None
        refq_ptr = None
        if ref_charges is not None:
            refq_arr = np.ascontiguousarray(ref_charges, dtype=np.float64)
            if refq_arr.shape != (len(at),):
                raise ValueError(
                    f"ref_charges must have shape ({len(at)},), got {refq_arr.shape}"
                )
            refq_ptr = refq_arr.ctypes.data_as(ctypes.c_void_p)

        bond_arr = None
        bond_ptr = None
        if bond_matrix is not None:
            bond_arr = np.ascontiguousarray(bond_matrix, dtype=np.int32)
            if bond_arr.shape != (len(at), len(at)):
                raise ValueError(
                    f"bond_matrix must have shape ({len(at)}, {len(at)}), "
                    f"got {bond_arr.shape}"
                )
            bond_ptr = bond_arr.ctypes.data_as(ctypes.c_void_p)

        version_int = _resolve_version(version)
        # The Fortran side signals a bad parameter file the same way it signals
        # any other setup failure -- a NULL handle. Check the two cases we can
        # name precisely here, so they do not surface as a generic init error.
        if parametrisation:
            if not os.path.isfile(parametrisation):
                raise FileNotFoundError(
                    f"parameter file not found: {parametrisation}"
                )
            if parametrisation.endswith(".toml") and not toml_available():
                raise RuntimeError(
                    f"cannot read {parametrisation}: libgfnff was built without "
                    "toml-f, so TOML parameter files are unsupported. Rebuild "
                    "with -DWITH_TOMLF=ON (CMake) or -Dtoml-f=enabled (meson)."
                )
        param_bytes = parametrisation.encode() if parametrisation else None

        use_ex = (
            fragments is not None
            or ref_charges is not None
            or accuracy is not None
            or version_int != VERSION_DEFAULT
            or param_bytes is not None
            or bond_ptr is not None
        )

        if use_ex:
            # extended path: carries the optional host-supplied bundle hints
            if lattice is not None:
                lat = np.ascontiguousarray(lattice, dtype=np.float64)
                if lat.shape != (3, 3):
                    raise ValueError(
                        f"lattice must have shape (3, 3), got {lat.shape}"
                    )
                lattice_ptr = lat.ctypes.data_as(ctypes.c_void_p)
            else:
                lattice_ptr = None
            self._handle = _lib.c_gfnff_calculator_init_ex(
                nat, at, xyz,
                ctypes.c_int(charge),
                ctypes.c_int(printlevel),
                solvent_bytes,
                lattice_ptr,
                ctypes.c_int(npbc),
                frag_ptr,
                refq_ptr,
                ctypes.c_double(-1.0 if accuracy is None else float(accuracy)),
                ctypes.c_int(version_int),
                param_bytes,
                bond_ptr,
            )
        elif lattice is not None:
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
            hint = ""
            if version_int == int(Version.mcgfnff2023) and npbc == 0:
                hint = (
                    " The mcgfnff2023 parametrisation is for molecular crystals "
                    "and is rejected for non-periodic systems."
                )
            raise RuntimeError(
                "GFN-FF initialization failed (c_gfnff_calculator_init returned NULL). "
                "Check printlevel > 0 for error details." + hint
            )
        self._alive = True
        self._nat = len(at)
        #: kept for introspection; the library has no getter for either
        self.version = Version(version_int) if version_int in _version_values.values() \
            else Version.default
        self.parametrisation = parametrisation

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

        if lattice is not None:
            lat = np.ascontiguousarray(lattice, dtype=np.float64)
            lattice_ptr = lat.ctypes.data_as(ctypes.c_void_p)
        else:
            lattice_ptr = None

        _lib.c_gfnff_calculator_singlepoint(
            ctypes.byref(self._handle),
            ctypes.c_int(nat),
            at,
            xyz,
            ctypes.byref(energy_out),
            gradient,
            sigma,
            lattice_ptr,
            ctypes.byref(iostat),
        )

        if iostat.value != 0:
            raise RuntimeError(
                f"GFN-FF singlepoint failed with iostat = {iostat.value}"
            )

        return float(energy_out.value), gradient, sigma

    # ------------------------------------------------------------------

    def hessian(self, numbers, positions_bohr, step=None):
        """Compute the Cartesian nuclear Hessian for the given geometry.

        Only available for non-periodic systems: the analytic terms have no
        periodic implementation yet, and finite-differencing them would
        silently ignore the periodic images.

        Parameters
        ----------
        numbers:
            Atomic numbers, shape ``(nat,)``.
        positions_bohr:
            Coordinates in Bohr, shape ``(nat, 3)``.
        step:
            Finite-difference step in Bohr for the terms not yet available in
            closed form.  ``None`` lets the library choose.

        Returns
        -------
        hessian : ndarray, shape ``(3 * nat, 3 * nat)``
            Second derivatives in Eh / Bohr**2.  Degree of freedom ``(c, a)``
            is at index ``3 * a + c``.  The matrix is symmetric.
        energy : float
            Total energy in Hartree.
        gradient : ndarray, shape ``(nat, 3)``
            Analytic gradient in Eh / Bohr.

        Raises
        ------
        NotImplementedError
            If the calculator was initialized for a periodic system.
        RuntimeError
            If the library reports any other failure.
        """
        self._check_alive()

        at = np.asarray(numbers, dtype=np.int32)
        xyz = np.ascontiguousarray(positions_bohr, dtype=np.float64)
        nat = len(at)

        hessian = np.zeros((3 * nat, 3 * nat), dtype=np.float64, order="C")
        gradient = np.zeros((nat, 3), dtype=np.float64, order="C")
        energy_out = ctypes.c_double(0.0)
        iostat = ctypes.c_int(0)

        _lib.c_gfnff_calculator_hessian(
            ctypes.byref(self._handle),
            ctypes.c_int(nat),
            at,
            xyz,
            hessian,
            ctypes.cast(ctypes.byref(energy_out), ctypes.c_void_p),
            gradient.ctypes.data_as(ctypes.c_void_p),
            ctypes.c_double(-1.0 if step is None else float(step)),
            ctypes.byref(iostat),
        )

        if iostat.value == PBC_UNSUPPORTED:
            raise NotImplementedError(
                "GFN-FF Hessian is not available for periodic systems"
            )
        if iostat.value != 0:
            raise RuntimeError(
                f"GFN-FF Hessian failed with iostat = {iostat.value}"
            )

        return hessian, float(energy_out.value), gradient

    # ------------------------------------------------------------------

    def charges(self):
        """Atomic partial charges from the last singlepoint.

        The EEQ charges are a by-product of the energy evaluation rather than
        a separate model, so they belong to the geometry most recently passed
        to :meth:`singlepoint` (or :meth:`hessian`, which runs one).

        Returns
        -------
        charges : ndarray, shape ``(nat,)``
            Partial charges in units of the elementary charge.  They sum to
            the total charge given at initialization.

        Raises
        ------
        RuntimeError
            If no singlepoint has been run on this calculator yet, or if the
            selected version has no charges at all (``harmonic2020`` returns
            before the EEQ solve).
        """
        self._check_alive()

        charges = np.zeros(self._nat, dtype=np.float64, order="C")
        iostat = ctypes.c_int(0)

        _lib.c_gfnff_calculator_charges(
            ctypes.byref(self._handle),
            ctypes.c_int(self._nat),
            charges,
            ctypes.byref(iostat),
        )

        if iostat.value == NO_CHARGES:
            if self.version == Version.harmonic2020:
                raise RuntimeError(
                    "the harmonic2020 parametrisation computes no charges: it "
                    "returns before the EEQ solve"
                )
            raise RuntimeError(
                "no GFN-FF charges available yet; run singlepoint() first"
            )
        if iostat.value != 0:
            raise RuntimeError(
                f"GFN-FF charge lookup failed with iostat = {iostat.value}"
            )

        return charges

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
