# This file is part of gfnff.
# SPDX-Identifier: LGPL-3.0-or-later
"""
Tests for the ASE GFNFF calculator (Angstrom / eV units).
"""

import numpy as np
import pytest

pytest.importorskip("ase", reason="ase not installed")

from ase.units import Bohr, Hartree

from gfnff import GFNFF

from .conftest import CAFFEINE_ENERGY_HARTREE, SIO2_ENERGY_HARTREE


# ---------------------------------------------------------------------------
# Energy and forces
# ---------------------------------------------------------------------------

class TestEnergy:
    def test_energy_caffeine(self, caffeine_ase):
        calc = GFNFF()
        caffeine_ase.calc = calc
        energy = caffeine_ase.get_potential_energy()
        expected = CAFFEINE_ENERGY_HARTREE * Hartree
        assert abs(energy - expected) < 1e-4  # eV tolerance

    def test_energy_pbc_sio2(self, sio2_ase):
        calc = GFNFF()
        sio2_ase.calc = calc
        energy = sio2_ase.get_potential_energy()
        expected = SIO2_ENERGY_HARTREE * Hartree
        assert abs(energy - expected) < 1e-3  # eV tolerance (6-digit ref)


class TestForces:
    def test_forces_shape(self, caffeine_ase):
        calc = GFNFF()
        caffeine_ase.calc = calc
        forces = caffeine_ase.get_forces()
        assert forces.shape == (len(caffeine_ase), 3)
        assert forces.dtype == np.float64

    def test_forces_units_finite_diff(self, caffeine_ase):
        """Forces should match negative finite-difference gradient (eV/Å)."""
        calc = GFNFF()
        caffeine_ase.calc = calc
        forces = caffeine_ase.get_forces()

        delta = 1e-3  # Angstrom
        atoms_fwd = caffeine_ase.copy()
        atoms_fwd.positions[0, 0] += delta
        atoms_fwd.calc = GFNFF()
        e_fwd = atoms_fwd.get_potential_energy()

        atoms_bwd = caffeine_ase.copy()
        atoms_bwd.positions[0, 0] -= delta
        atoms_bwd.calc = GFNFF()
        e_bwd = atoms_bwd.get_potential_energy()

        fd_force_x = -(e_fwd - e_bwd) / (2 * delta)
        assert abs(fd_force_x - forces[0, 0]) < 1e-3  # eV/Å

    def test_forces_pbc_sio2(self, sio2_ase):
        calc = GFNFF()
        sio2_ase.calc = calc
        forces = sio2_ase.get_forces()
        assert forces.shape == (len(sio2_ase), 3)


# ---------------------------------------------------------------------------
# Stress tensor
# ---------------------------------------------------------------------------

class TestStress:
    def test_stress_shape(self, caffeine_ase):
        """Stress is always a 6-element Voigt array [xx,yy,zz,yz,xz,xy]."""
        calc = GFNFF()
        caffeine_ase.calc = calc
        stress = caffeine_ase.get_stress()
        assert stress.shape == (6,)
        assert stress.dtype == np.float64

    def test_stress_zero_for_nonpbc(self, caffeine_ase):
        """Non-periodic systems have zero stress."""
        calc = GFNFF()
        caffeine_ase.calc = calc
        stress = caffeine_ase.get_stress()
        assert np.allclose(stress, 0.0)

    def test_stress_nonzero_for_pbc(self, sio2_ase):
        """Periodic SiO2 must yield a non-zero stress tensor."""
        calc = GFNFF()
        sio2_ase.calc = calc
        stress = sio2_ase.get_stress()
        assert stress.shape == (6,)
        assert not np.allclose(stress, 0.0)


# ---------------------------------------------------------------------------
# Calculator state management
# ---------------------------------------------------------------------------

class TestStateManagement:
    def test_topology_reinit_on_number_change(self, caffeine_ase):
        """Changing atomic numbers should trigger re-initialization."""
        calc = GFNFF()
        caffeine_ase.calc = calc
        caffeine_ase.get_potential_energy()
        old_gfnff = calc._gfnff

        # Modify numbers (e.g. change last H to D — same mass, same number
        # so tweak a carbon to nitrogen instead)
        caffeine_ase.numbers[0] = 7  # C → N
        caffeine_ase.get_potential_energy()
        assert calc._gfnff is not old_gfnff, "topology should have been rebuilt"

    def test_topology_reused_for_position_change(self, caffeine_ase):
        """Topology should NOT be rebuilt for a pure geometry update."""
        calc = GFNFF()
        caffeine_ase.calc = calc
        caffeine_ase.get_potential_energy()
        old_gfnff = calc._gfnff

        caffeine_ase.positions[0, 0] += 0.01  # tiny displacement
        caffeine_ase.get_potential_energy()
        assert calc._gfnff is old_gfnff, "topology should have been reused"

    def test_charge_from_atoms_info(self, caffeine_ase):
        """Charge in atoms.info overrides the calculator parameter."""
        calc = GFNFF(charge=0)
        caffeine_ase.info["charge"] = 0
        caffeine_ase.calc = calc
        e_neutral = caffeine_ase.get_potential_energy()

        caffeine_ase.info["charge"] = 1
        caffeine_ase2 = caffeine_ase.copy()
        caffeine_ase2.info["charge"] = 1
        caffeine_ase2.calc = GFNFF(charge=0)
        e_charged = caffeine_ase2.get_potential_energy()

        assert abs(e_neutral - e_charged) > 0.01  # eV — charge changes energy


# ---------------------------------------------------------------------------
# Solvation
# ---------------------------------------------------------------------------

class TestSolvation:
    def test_solvent_changes_energy(self, caffeine_ase):
        calc_vac = GFNFF(solvent="")
        caffeine_ase.calc = calc_vac
        e_vac = caffeine_ase.get_potential_energy()

        caffeine_ase2 = caffeine_ase.copy()
        calc_sol = GFNFF(solvent="h2o")
        caffeine_ase2.calc = calc_sol
        e_sol = caffeine_ase2.get_potential_energy()

        assert e_sol != e_vac, "solvation should change the energy"
