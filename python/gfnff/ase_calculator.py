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

from .calculator import GFNFFCalculator


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
    """

    implemented_properties = ["energy", "forces", "stress"]

    default_parameters = {
        "charge": 0,
        "solvent": "",
        "printlevel": 0,
    }

    def __init__(self, charge=0, solvent="", printlevel=0, **kwargs):
        super().__init__(**kwargs)
        self.parameters.update(
            charge=charge,
            solvent=solvent,
            printlevel=printlevel,
        )
        self._gfnff: GFNFFCalculator | None = None
        self._last_numbers: np.ndarray | None = None
        self._last_pbc: np.ndarray | None = None

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _needs_reinit(self, atoms, system_changes) -> bool:
        """Topology must be rebuilt when atom types or PBC/cell change."""
        if self._gfnff is None:
            return True
        if not np.array_equal(atoms.numbers, self._last_numbers):
            return True
        if "cell" in system_changes or "pbc" in system_changes:
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

        return GFNFFCalculator(
            numbers,
            pos_bohr,
            charge=int(charge),
            printlevel=self.parameters.printlevel,
            solvent=self.parameters.solvent,
            lattice=lattice_bohr,
            npbc=npbc,
        )

    # ------------------------------------------------------------------
    # Calculator interface
    # ------------------------------------------------------------------

    def calculate(self, atoms=None, properties=None, system_changes=all_changes):
        if properties is None:
            properties = self.implemented_properties

        Calculator.calculate(self, atoms, properties, system_changes)

        atoms = self.atoms

        if self._needs_reinit(atoms, system_changes):
            if self._gfnff is not None:
                self._gfnff.deallocate()
            self._gfnff = self._make_gfnff(atoms)
            self._last_numbers = atoms.numbers.copy()
            self._last_pbc = atoms.pbc.copy()

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
