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


# ---------------------------------------------------------------------------
# Partial charges
# ---------------------------------------------------------------------------

class TestCharges:
    def test_get_charges(self, caffeine_ase):
        caffeine_ase.calc = GFNFF()
        q = caffeine_ase.get_charges()
        assert q.shape == (len(caffeine_ase),)
        assert q.sum() == pytest.approx(0.0, abs=1e-8)
        numbers = caffeine_ase.numbers
        assert (q[numbers == 8] < 0).all()
        assert (q[numbers == 1] > 0).all()

    def test_charges_are_a_declared_property(self):
        assert "charges" in GFNFF.implemented_properties

    def test_charges_track_geometry(self, caffeine_ase):
        caffeine_ase.calc = GFNFF()
        first = caffeine_ase.get_charges().copy()
        caffeine_ase.positions[0, 0] += 0.15  # Angstrom
        second = caffeine_ase.get_charges()
        assert np.abs(first - second).max() > 1e-6

    def test_charges_are_unscaled(self, caffeine_ase):
        """Charges are in e and must not pick up any eV/Angstrom conversion."""
        from gfnff import GFNFFCalculator

        caffeine_ase.calc = GFNFF()
        ase_q = caffeine_ase.get_charges()

        with GFNFFCalculator(
            np.asarray(caffeine_ase.numbers, dtype=np.int32),
            np.ascontiguousarray(caffeine_ase.positions / Bohr),
        ) as raw:
            raw.singlepoint(
                np.asarray(caffeine_ase.numbers, dtype=np.int32),
                np.ascontiguousarray(caffeine_ase.positions / Bohr),
            )
            assert np.allclose(ase_q, raw.charges(), atol=0.0, rtol=0.0)

    def test_charges_pbc(self, sio2_ase):
        sio2_ase.calc = GFNFF()
        q = sio2_ase.get_charges()
        assert q.sum() == pytest.approx(0.0, abs=1e-8)


# ---------------------------------------------------------------------------
# Parametrisation selection
# ---------------------------------------------------------------------------

class TestVersion:
    def test_harmonic_changes_energy(self, caffeine_ase):
        caffeine_ase.calc = GFNFF()
        default = caffeine_ase.get_potential_energy()

        other = caffeine_ase.copy()
        other.calc = GFNFF(version="harmonic2020")
        harmonic = other.get_potential_energy()

        assert abs(harmonic - default) > 1e-6

    def test_harmonic_omits_charges(self, caffeine_ase):
        """No EEQ solve runs, so the property is absent rather than zeroed."""
        from ase.calculators.calculator import PropertyNotImplementedError

        caffeine_ase.calc = GFNFF(version="harmonic2020")
        assert np.isfinite(caffeine_ase.get_potential_energy())
        with pytest.raises((PropertyNotImplementedError, RuntimeError)):
            caffeine_ase.get_charges()

    def test_version_change_rebuilds_topology(self, caffeine_ase):
        """set() must not leave the old force field in place."""
        calc = GFNFF()
        caffeine_ase.calc = calc
        default = caffeine_ase.get_potential_energy()

        calc.set(version="harmonic2020")
        harmonic = caffeine_ase.get_potential_energy()

        assert abs(harmonic - default) > 1e-6

    def test_conformer_matches_default(self, caffeine_ase):
        """Near equilibrium the two share the same bond well; see test_lib."""
        caffeine_ase.calc = GFNFF()
        default = caffeine_ase.get_potential_energy()

        other = caffeine_ase.copy()
        other.calc = GFNFF(version="conformer2020")
        assert abs(other.get_potential_energy() - default) < 1e-10

    def test_conformer_bond_cannot_dissociate(self, caffeine_ase):
        """One calculator, then move an atom: the topology stays put.

        This is the path an optimiser or an MD run takes, and the only one
        where the bond term is still asked about a bond that has been pulled
        apart -- rebuilding the calculator would simply drop the bond.
        """
        direction = caffeine_ase.positions[14] - caffeine_ase.positions[0]
        direction /= np.linalg.norm(direction)

        def curve(**kwargs):
            atoms = caffeine_ase.copy()
            atoms.calc = GFNFF(**kwargs)
            origin = atoms.positions[14].copy()
            out = []
            for k in range(4):
                atoms.positions[14] = origin + direction * k
                out.append(atoms.get_potential_energy())
            return out

        ref = curve()
        con = curve(version="conformer2020")
        assert (ref[3] - ref[2]) / (ref[2] - ref[1]) < 0.1
        assert 0.9 < (con[3] - con[2]) / (con[2] - con[1]) < 1.1
        assert con[3] > ref[3] + 1.0

    def test_mcgfnff_for_crystal(self, sio2_ase):
        sio2_ase.calc = GFNFF(version="mcgfnff2023")
        assert np.isfinite(sio2_ase.get_potential_energy())

    def test_parameter_change_discards_cached_results(self, caffeine_ase):
        """ASE's base Calculator caches across set(); this one must not."""
        calc = GFNFF()
        caffeine_ase.calc = calc
        vacuum = caffeine_ase.get_potential_energy()

        calc.set(solvent="h2o")
        solvated = caffeine_ase.get_potential_energy()

        assert solvated != vacuum


# ---------------------------------------------------------------------------
# Hessian
# ---------------------------------------------------------------------------

class TestHessian:
    def test_shape_and_symmetry(self, caffeine_ase):
        caffeine_ase.calc = GFNFF()
        h = caffeine_ase.calc.get_hessian(caffeine_ase)
        n3 = 3 * len(caffeine_ase)
        assert h.shape == (n3, n3)
        assert np.abs(h - h.T).max() == 0.0

    def test_units_against_finite_differenced_forces(self, caffeine_ase):
        """The whole point of the ASE layer is the unit conversion.

        Differencing ASE forces (eV/A) over an ASE displacement (A) gives
        eV/A**2 directly, so this fails if the Hartree or Bohr factor is
        missing, doubled, or inverted -- none of which a shape or symmetry
        check would catch.
        """
        caffeine_ase.calc = GFNFF()
        h = caffeine_ase.calc.get_hessian(caffeine_ase)

        delta = 1.0e-3  # Angstrom
        shifted = caffeine_ase.copy()
        shifted.calc = GFNFF()
        shifted.positions[0, 0] += delta
        fplus = shifted.get_forces()
        shifted.positions[0, 0] -= 2.0 * delta
        fminus = shifted.get_forces()

        # forces are -gradient, so the second derivative picks up a sign
        fd = (-(fplus - fminus) / (2.0 * delta)).reshape(-1)
        assert np.abs(fd - h[:, 0]).max() < 1.0e-3

    def test_not_an_implemented_property(self, caffeine_ase):
        """It must not ride along on every calculate() call."""
        assert "hessian" not in GFNFF.implemented_properties

    def test_rejects_periodic(self, sio2_ase):
        sio2_ase.calc = GFNFF()
        with pytest.raises(NotImplementedError, match="periodic"):
            sio2_ase.calc.get_hessian(sio2_ase)

    def test_uses_attached_atoms_by_default(self, caffeine_ase):
        caffeine_ase.calc = GFNFF()
        caffeine_ase.get_potential_energy()  # attaches atoms
        explicit = caffeine_ase.calc.get_hessian(caffeine_ase)
        implicit = caffeine_ase.calc.get_hessian()
        assert np.array_equal(explicit, implicit)

    def test_cache_is_geometry_specific(self, caffeine_ase):
        calc = GFNFF()
        caffeine_ase.calc = calc
        first = calc.get_hessian(caffeine_ase)
        moved = caffeine_ase.copy()
        moved.positions[0, 0] += 0.2
        second = calc.get_hessian(moved)
        assert np.abs(first - second).max() > 1e-6

    def test_cache_does_not_leak_a_mutable_reference(self, caffeine_ase):
        calc = GFNFF()
        caffeine_ase.calc = calc
        first = calc.get_hessian(caffeine_ase)
        first[0, 0] = 12345.0
        second = calc.get_hessian(caffeine_ase)
        assert second[0, 0] != 12345.0

    def test_step_argument_is_accepted(self, caffeine_ase):
        calc = GFNFF()
        caffeine_ase.calc = calc
        default = calc.get_hessian(caffeine_ase)
        calc.reset()
        tighter = calc.get_hessian(caffeine_ase, step=1.0e-3 * Bohr)
        assert np.abs(default - tighter).max() < 1.0e-2

    def test_vibrations_data(self, caffeine_ase):
        """Frequencies through ASE's own machinery, no 6N displacements."""
        pytest.importorskip("ase.vibrations")
        caffeine_ase.calc = GFNFF()
        vib = caffeine_ase.calc.get_vibrations(caffeine_ase)
        energies = vib.get_energies()
        assert len(energies) == 3 * len(caffeine_ase)
        # an equilibrium-ish structure: at most the 6 translations/rotations
        # may come out imaginary
        assert (energies.imag > 0).sum() <= 6


class TestBondMatrix:
    """Host-supplied molecular graph through the ASE calculator."""

    @staticmethod
    def _graph(atoms, cutoff=1.7):
        """Adjacency from distances in Angstrom; exact for caffeine."""
        d = atoms.get_all_distances()
        n = len(atoms)
        bm = ((d < cutoff) & (d > 0.0)).astype(np.int32)
        assert bm.sum() // 2 > 0
        return bm

    def test_graph_reproduces_perception(self, caffeine_ase):
        caffeine_ase.calc = GFNFF()
        default = caffeine_ase.get_potential_energy()

        other = caffeine_ase.copy()
        other.calc = GFNFF(bond_matrix=self._graph(caffeine_ase))
        assert abs(other.get_potential_energy() - default) < 1e-10

    def test_graph_via_atoms_info(self, caffeine_ase):
        caffeine_ase.calc = GFNFF()
        default = caffeine_ase.get_potential_energy()

        other = caffeine_ase.copy()
        other.info["bond_matrix"] = self._graph(caffeine_ase)
        other.calc = GFNFF()
        assert abs(other.get_potential_energy() - default) < 1e-10

    def test_different_graph_is_a_different_force_field(self, caffeine_ase):
        """Dropping a bond from the graph must change the energy."""
        bm = self._graph(caffeine_ase)
        stripped = bm.copy()
        stripped[0, 14] = stripped[14, 0] = 0   # drop one C-H

        full = caffeine_ase.copy()
        full.calc = GFNFF(bond_matrix=bm)
        less = caffeine_ase.copy()
        less.calc = GFNFF(bond_matrix=stripped)
        assert abs(less.get_potential_energy() - full.get_potential_energy()) > 1e-6

    def test_graph_is_part_of_the_rebuild_key(self, caffeine_ase):
        """The graph defines the topology, so it has to force a rebuild.

        Like ``fragments`` and ``ref_charges``, it is a constructor argument
        rather than a ``set()``-able parameter, so this checks the mechanism
        that decides a rebuild rather than going through ``set()``.
        """
        bm = self._graph(caffeine_ase)
        calc = GFNFF(bond_matrix=bm)
        caffeine_ase.calc = calc
        before = calc._setup_key(caffeine_ase)

        stripped = bm.copy()
        stripped[0, 14] = stripped[14, 0] = 0
        caffeine_ase.info["bond_matrix"] = stripped
        assert calc._setup_key(caffeine_ase) != before

    def test_soup_to_structure(self, caffeine_ase):
        """The application: a cloud of atoms plus a graph becomes a molecule.

        Bond perception on the scrambled coordinates would find a different
        molecule entirely, so this only works because the graph replaces it.
        """
        from ase.optimize import BFGS

        ref = caffeine_ase.copy()
        ref.calc = GFNFF()
        e_ref = ref.get_potential_energy()
        bm = self._graph(caffeine_ase)

        soup = caffeine_ase.copy()
        rng = np.random.default_rng(20240803)
        soup.positions = rng.uniform(-6.0, 6.0, size=(len(soup), 3))

        soup.calc = GFNFF(version="harmonic2020", bond_matrix=bm)
        BFGS(soup, logfile=None).run(fmax=1e-3, steps=3000)

        # every declared bond must have come out at a chemically sane length
        bonded = np.array([
            soup.get_distance(i, j)
            for i in range(len(soup)) for j in range(i) if bm[i, j]
        ])
        assert bonded.max() < 2.0, f"longest recovered bond {bonded.max():.2f} A"
        assert bonded.min() > 0.7, f"shortest recovered bond {bonded.min():.2f} A"

        # and regular GFN-FF must then relax it to a real minimum
        soup.calc = GFNFF()
        BFGS(soup, logfile=None).run(fmax=1e-3, steps=3000)
        assert soup.get_potential_energy() < e_ref + 1.0
