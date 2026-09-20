# This file is part of gfnff.
# SPDX-Identifier: LGPL-3.0-or-later
"""
ASE Calculator interface for GFN-FF.

Unit conventions
----------------
ASE uses Angstrom and eV throughout.  The conversions applied here are:

  positions   : Angstrom  →  Bohr       (divide by ase.units.Bohr)
  lattice     : Angstrom  →  Bohr       (divide by ase.units.Bohr)
  energy      : Hartree   →  eV         (multiply by ase.units.Hartree)
  forces      : -(Eh/Bohr) → eV/Ang    (multiply by ase.units.Hartree / ase.units.Bohr)
  stress      : Eh/Bohr³  → eV/Ang³    (sigma / volume * Hartree / Bohr**3)
                returned in Voigt order [xx, yy, zz, yz, xz, xy]

For non-periodic systems sigma is zero, so stress is reported as a zero
six-vector.  For periodic systems the stress is sigma / cell_volume.
"""

import numpy as np

from ase.calculators.calculator import Calculator, all_changes
from ase.units import Bohr, Hartree

from .calculator import GFNFFCalculator, Version


def _graph_fingerprint(bond_matrix):
    """A hashable stand-in for a bond matrix, for change detection.

    ``_setup_key`` compares two dicts with ``!=``; an array in there would
    raise on the elementwise comparison, so the graph is reduced to bytes.
    """
    if bond_matrix is None:
        return None
    arr = np.ascontiguousarray(bond_matrix, dtype=np.int32)
    return (arr.shape, arr.tobytes())


class GFNFF(Calculator):
    """ASE calculator for the GFN-FF force field.

    Parameters
    ----------
    charge : int, optional
        Total molecular charge.  May also be provided via
        ``atoms.info["charge"]``, which takes precedence.
    solvent : str, optional
        Implicit solvent name (e.g. ``"h2o"``, ``"acetone"``).
        Empty string (default) runs in vacuum.
    printlevel : int, optional
        Fortran output verbosity (0 = silent).
    fragments : array_like of int, optional
        Per-atom fragment index, shape ``(nat,)``.  When given, GFN-FF forms no
        bonds between atoms of differing fragments (host-defined fragmentation
        instead of automatic detection).  May also be provided via
        ``atoms.info["fragments"]``, which takes precedence.
    ref_charges : array_like of float, optional
        Per-atom reference charges, shape ``(nat,)``, summed per fragment to set
        the per-fragment EEQ charge constraint.  No charge model is invoked.
        May also be provided via ``atoms.info["ref_charges"]`` (takes precedence).
    version : Version or int or str, optional
        Parametrisation version, e.g. ``"mcgfnff2023"``.  ``None`` uses the
        library default.  ``mcgfnff2023`` requires a periodic system.
    parametrisation : str, optional
        Path to a parameter file.  A ``.toml`` path is applied as an overlay on
        the internal set for ``version``, so it need only name the keys it
        changes.
    bond_matrix : array_like of int, optional
        Molecular graph, shape ``(nat, nat)`` of integer bond orders.  Replaces
        GFN-FF's distance-based bond perception, so the connectivity need not
        be supported by the coordinates.  Combined with ``version="harmonic2020"``
        this turns an unstructured cloud of atoms plus a graph into a 3D
        structure.  May also be provided via ``atoms.info["bond_matrix"]``,
        which takes precedence.

    Notes
    -----
    ``atoms.get_charges()`` returns the EEQ partial charges from the current
    geometry.  They are a by-product of the energy evaluation, so asking for
    them triggers a calculation exactly as energy or forces would.
    """

    implemented_properties = ["energy", "forces", "stress", "charges"]

    # Every parameter this calculator takes -- charge, solvent, version,
    # parametrisation, accuracy -- changes the numbers. ASE's base Calculator
    # leaves this False, which would let set() swap the force field out and
    # still hand back the cached energy from the old one.
    discard_results_on_any_change = True

    default_parameters = {
        "charge": 0,
        "solvent": "",
        "printlevel": 0,
        "accuracy": None,
        "version": None,
        "parametrisation": None,
    }

    def __init__(self, charge=0, solvent="", printlevel=0,
                 fragments=None, ref_charges=None, accuracy=None,
                 version=None, parametrisation=None, bond_matrix=None,
                 **kwargs):
        super().__init__(**kwargs)
        self.parameters.update(
            charge=charge,
            solvent=solvent,
            printlevel=printlevel,
            accuracy=accuracy,
            version=version,
            parametrisation=parametrisation,
        )
        # per-atom arrays are kept as attributes (not in self.parameters, which
        # is meant for simple, serializable values used in change detection)
        self._fragments = fragments
        self._ref_charges = ref_charges
        self._bond_matrix = bond_matrix
        self._gfnff: GFNFFCalculator | None = None
        self._last_numbers: np.ndarray | None = None
        self._last_pbc: np.ndarray | None = None
        self._last_setup: dict | None = None
        self._hessian_cache: tuple[np.ndarray, np.ndarray] | None = None

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _setup_key(self, atoms) -> dict:
        """Everything the Fortran calculator is built from, except geometry.

        Charge, solvent, accuracy, version and the parameter file are all fixed
        when the calculator is constructed, so changing any of them through
        ``set()`` has to rebuild it.  Discarding the cached results is not
        enough: the next singlepoint would otherwise run on the old setup and
        quietly return the old answer under the new parameters.
        """
        return {
            "charge": int(atoms.info.get("charge", self.parameters.charge)),
            "solvent": self.parameters.solvent,
            "accuracy": self.parameters.accuracy,
            "version": self.parameters.version,
            "parametrisation": self.parameters.parametrisation,
            # the graph defines the topology, so swapping it is as much a
            # rebuild as swapping the parametrisation. Hashed rather than
            # carried, since the key is compared with != on a plain dict.
            "bond_matrix": _graph_fingerprint(
                atoms.info.get("bond_matrix", self._bond_matrix)
            ),
        }

    def _needs_reinit(self, atoms, system_changes) -> bool:
        """Topology must be rebuilt when atom types or PBC flags change.
        Cell-only changes are handled by passing the updated lattice to
        singlepoint, which updates the cell in-place without redoing the
        full topology setup."""
        if self._gfnff is None:
            return True
        if not np.array_equal(atoms.numbers, self._last_numbers):
            return True
        if "pbc" in system_changes and not np.array_equal(atoms.pbc, self._last_pbc):
            return True
        if self._setup_key(atoms) != self._last_setup:
            return True
        return False

    def _make_gfnff(self, atoms) -> GFNFFCalculator:
        numbers = np.asarray(atoms.numbers, dtype=np.int32)
        pos_bohr = np.ascontiguousarray(atoms.positions / Bohr, dtype=np.float64)
        charge = atoms.info.get("charge", self.parameters.charge)

        npbc = int(np.sum(atoms.pbc))
        if npbc > 0:
            lattice_bohr = np.ascontiguousarray(
                atoms.cell[:] / Bohr, dtype=np.float64
            )
        else:
            lattice_bohr = None

        # optional host-supplied hints; atoms.info overrides constructor values
        fragments = atoms.info.get("fragments", self._fragments)
        ref_charges = atoms.info.get("ref_charges", self._ref_charges)

        return GFNFFCalculator(
            numbers,
            pos_bohr,
            charge=int(charge),
            printlevel=self.parameters.printlevel,
            solvent=self.parameters.solvent,
            lattice=lattice_bohr,
            npbc=npbc,
            fragments=fragments,
            ref_charges=ref_charges,
            accuracy=self.parameters.accuracy,
            version=self.parameters.version,
            parametrisation=self.parameters.parametrisation,
            bond_matrix=atoms.info.get("bond_matrix", self._bond_matrix),
        )

    def _ensure_gfnff(self, atoms, system_changes):
        """Build the Fortran calculator, or keep the existing one."""
        if self._needs_reinit(atoms, system_changes):
            if self._gfnff is not None:
                self._gfnff.deallocate()
            self._gfnff = self._make_gfnff(atoms)
            self._last_numbers = atoms.numbers.copy()
            self._last_pbc = atoms.pbc.copy()
            self._last_setup = self._setup_key(atoms)

    # ------------------------------------------------------------------
    # Calculator interface
    # ------------------------------------------------------------------

    def calculate(self, atoms=None, properties=None, system_changes=all_changes):
        if properties is None:
            properties = self.implemented_properties

        Calculator.calculate(self, atoms, properties, system_changes)

        atoms = self.atoms

        self._ensure_gfnff(atoms, system_changes)

        # Prepare geometry for this step
        numbers = np.asarray(atoms.numbers, dtype=np.int32)
        pos_bohr = np.ascontiguousarray(atoms.positions / Bohr, dtype=np.float64)

        npbc = int(np.sum(atoms.pbc))
        if npbc > 0:
            lattice_bohr = np.ascontiguousarray(
                atoms.cell[:] / Bohr, dtype=np.float64
            )
        else:
            lattice_bohr = None

        energy_ha, gradient, sigma = self._gfnff.singlepoint(
            numbers, pos_bohr, lattice=lattice_bohr
        )

        # Convert to ASE units
        self.results["energy"] = energy_ha * Hartree
        # forces = -gradient; convert Eh/Bohr → eV/Ang
        self.results["forces"] = -gradient * (Hartree / Bohr)
        # stress: sigma (Eh) / volume (Bohr³) → eV/Å³, Voigt order [xx,yy,zz,yz,xz,xy]
        if npbc > 0:
            volume_bohr3 = np.linalg.det(lattice_bohr)
            stress_3x3 = sigma / volume_bohr3 * (Hartree / Bohr**3)
        else:
            stress_3x3 = np.zeros((3, 3), dtype=np.float64)
        self.results["stress"] = stress_3x3[
            [0, 1, 2, 1, 0, 0], [0, 1, 2, 2, 2, 1]
        ]
        # EEQ partial charges, in e; they fall out of the singlepoint above and
        # need no unit conversion. harmonic2020 runs no EEQ solve, so the
        # property is simply left absent there rather than reported as zeros.
        try:
            self.results["charges"] = self._gfnff.charges()
        except RuntimeError:
            self.results.pop("charges", None)

    # ------------------------------------------------------------------
    # Hessian
    # ------------------------------------------------------------------

    def get_hessian(self, atoms=None, step=None):
        """Cartesian nuclear Hessian in eV / Angstrom**2.

        Deliberately *not* one of ``implemented_properties``: ASE would then
        compute it on every ``calculate`` call, and an MD step does not need a
        Hessian.  Ask for it explicitly instead.

        Parameters
        ----------
        atoms:
            Structure to evaluate.  Defaults to the one already attached.
        step:
            Finite-difference step in Angstrom for the terms not yet available
            in closed form.  ``None`` lets the library choose.

        Returns
        -------
        hessian : ndarray, shape ``(3 * nat, 3 * nat)``
            Second derivatives in eV / Angstrom**2.  Degree of freedom
            ``(c, a)`` is at index ``3 * a + c``.  The matrix is symmetric.

        Raises
        ------
        NotImplementedError
            For periodic systems: the analytic terms have no periodic
            implementation, and finite-differencing them would silently ignore
            the images.
        """
        # self.atoms is deliberately left alone: it is paired with self.results,
        # and repointing it at another structure would leave the cached energy
        # and forces describing a system the calculator no longer claims.
        atoms = self.atoms if atoms is None else atoms
        if atoms is None:
            raise ValueError("get_hessian() needs an Atoms object")

        if atoms.pbc.any():
            raise NotImplementedError(
                "GFN-FF Hessian is not available for periodic systems"
            )

        numbers = np.asarray(atoms.numbers, dtype=np.int32)
        positions = np.ascontiguousarray(atoms.positions, dtype=np.float64)
        # Repeat calls on an unchanged geometry are common -- frequencies and
        # normal modes are usually wanted together -- and a Hessian is far too
        # expensive to recompute for them. Numbers are part of the key: same
        # coordinates with different elements is a different Hessian.
        if self._hessian_cache is not None:
            key_numbers, key_positions, cached = self._hessian_cache
            if np.array_equal(key_numbers, numbers) and np.array_equal(
                key_positions, positions
            ):
                return cached.copy()

        self._ensure_gfnff(atoms, all_changes)

        pos_bohr = np.ascontiguousarray(positions / Bohr, dtype=np.float64)
        step_bohr = None if step is None else float(step) / Bohr

        hessian_au, _, _ = self._gfnff.hessian(numbers, pos_bohr, step=step_bohr)

        # Eh/Bohr**2 -> eV/Angstrom**2
        hessian = hessian_au * (Hartree / Bohr**2)
        self._hessian_cache = (numbers.copy(), positions.copy(), hessian)
        return hessian.copy()

    def get_vibrations(self, atoms=None, step=None):
        """Hessian wrapped as an ASE ``VibrationsData``.

        Gives frequencies, normal modes and thermochemistry through ASE's own
        machinery, without the finite-difference sweep ``ase.vibrations`` would
        otherwise run over 6N displaced singlepoints.
        """
        from ase.vibrations import VibrationsData

        atoms = self.atoms if atoms is None else atoms
        hessian = self.get_hessian(atoms=atoms, step=step)
        nat = len(atoms)
        # VibrationsData wants (nat, 3, nat, 3), not (3*nat, 3*nat)
        return VibrationsData(atoms, hessian.reshape(nat, 3, nat, 3))

    def reset(self):
        super().reset()
        self._hessian_cache = None
