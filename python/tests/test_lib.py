# This file is part of gfnff.
# SPDX-Identifier: LGPL-3.0-or-later
"""
Unit tests for the low-level GFNFFCalculator (Bohr / Hartree units).
These mirror the test cases in test/main.c.
"""

import numpy as np
import pytest

from gfnff import GFNFFCalculator, Version
from gfnff._lib import toml_available

from .conftest import CAFFEINE_ENERGY_HARTREE, SIO2_ENERGY_HARTREE


class TestInit:
    def test_init_caffeine(self, caffeine):
        numbers, xyz = caffeine
        with GFNFFCalculator(numbers, xyz) as calc:
            assert calc._handle.ptr is not None

    def test_init_pbc_sio2(self, sio2):
        numbers, xyz, lattice, npbc = sio2
        with GFNFFCalculator(numbers, xyz, lattice=lattice, npbc=npbc) as calc:
            assert calc._handle.ptr is not None

    def test_bad_positions_shape(self, caffeine):
        numbers, xyz = caffeine
        with pytest.raises(ValueError, match="shape"):
            GFNFFCalculator(numbers, xyz.flatten())

    def test_mismatched_length(self, caffeine):
        numbers, xyz = caffeine
        with pytest.raises(ValueError, match="same length"):
            GFNFFCalculator(numbers[:5], xyz)


class TestSinglepoint:
    def test_energy_caffeine(self, caffeine):
        numbers, xyz = caffeine
        with GFNFFCalculator(numbers, xyz) as calc:
            energy, gradient, sigma = calc.singlepoint(numbers, xyz)

        assert abs(energy - CAFFEINE_ENERGY_HARTREE) < 1e-6
        assert gradient.shape == (len(numbers), 3)
        assert gradient.dtype == np.float64
        assert sigma.shape == (3, 3)

    def test_sigma_zero_for_nonpbc(self, caffeine):
        """Stress tensor must be zero for non-periodic systems."""
        numbers, xyz = caffeine
        with GFNFFCalculator(numbers, xyz) as calc:
            _, _, sigma = calc.singlepoint(numbers, xyz)
        assert np.allclose(sigma, 0.0)

    def test_energy_pbc_sio2(self, sio2):
        numbers, xyz, lattice, npbc = sio2
        with GFNFFCalculator(numbers, xyz, lattice=lattice, npbc=npbc) as calc:
            energy, gradient, sigma = calc.singlepoint(numbers, xyz, lattice=lattice)

        assert abs(energy - SIO2_ENERGY_HARTREE) < 1e-4
        assert gradient.shape == (len(numbers), 3)
        assert sigma.shape == (3, 3)
        assert sigma.dtype == np.float64

    def test_sigma_nonzero_for_pbc(self, sio2):
        """Stress tensor must be non-zero for periodic systems."""
        numbers, xyz, lattice, npbc = sio2
        with GFNFFCalculator(numbers, xyz, lattice=lattice, npbc=npbc) as calc:
            _, _, sigma = calc.singlepoint(numbers, xyz, lattice=lattice)
        assert not np.allclose(sigma, 0.0)

    def test_gradient_finite_difference(self, caffeine):
        """First atom x-gradient roughly matches finite difference."""
        numbers, xyz = caffeine
        delta = 1e-4  # Bohr
        with GFNFFCalculator(numbers, xyz) as calc:
            _, grad, _ = calc.singlepoint(numbers, xyz)
            xyz_fwd = xyz.copy(); xyz_fwd[0, 0] += delta
            e_fwd, _, _ = calc.singlepoint(numbers, xyz_fwd)
            xyz_bwd = xyz.copy(); xyz_bwd[0, 0] -= delta
            e_bwd, _, _ = calc.singlepoint(numbers, xyz_bwd)

        fd = (e_fwd - e_bwd) / (2 * delta)
        assert abs(fd - grad[0, 0]) < 1e-4


class TestLifetime:
    def test_context_manager(self, caffeine):
        numbers, xyz = caffeine
        with GFNFFCalculator(numbers, xyz) as calc:
            assert calc._alive
        assert not calc._alive

    def test_double_deallocate(self, caffeine):
        numbers, xyz = caffeine
        calc = GFNFFCalculator(numbers, xyz)
        calc.deallocate()
        calc.deallocate()  # must not raise or crash

    def test_use_after_deallocate(self, caffeine):
        numbers, xyz = caffeine
        calc = GFNFFCalculator(numbers, xyz)
        calc.deallocate()
        with pytest.raises(RuntimeError, match="deallocated"):
            calc.singlepoint(numbers, xyz)

    def test_print_results(self, caffeine, capsys):
        numbers, xyz = caffeine
        with GFNFFCalculator(numbers, xyz) as calc:
            calc.singlepoint(numbers, xyz)
            calc.print_results()
        # Fortran writes directly to stdout unit 6; captured output may be empty
        # in some environments, so we just verify no exception is raised.


class TestHessian:
    """The Cartesian nuclear Hessian, exposed through c_gfnff_calculator_hessian."""

    def test_shape_and_symmetry(self, caffeine):
        numbers, positions = caffeine
        calc = GFNFFCalculator(numbers, positions)
        hess, energy, gradient = calc.hessian(numbers, positions)

        n3 = 3 * len(numbers)
        assert hess.shape == (n3, n3)
        assert gradient.shape == (len(numbers), 3)
        assert energy == pytest.approx(CAFFEINE_ENERGY_HARTREE, abs=1e-6)
        # the library symmetrises, so any asymmetry means the row-major
        # handoff from Fortran is wrong
        assert np.abs(hess - hess.T).max() == 0.0

    def test_matches_finite_difference(self, caffeine):
        """One column against a central difference of the API's own gradient.

        Checks the whole path, not just the shape: a transposed or misindexed
        buffer would still be symmetric and still have the right norm.
        """
        numbers, positions = caffeine
        calc = GFNFFCalculator(numbers, positions)
        hess, _, _ = calc.hessian(numbers, positions)

        delta = 5.0e-4
        shifted = positions.copy()
        shifted[0, 0] += delta
        _, gplus, _ = calc.singlepoint(numbers, shifted)
        shifted[0, 0] -= 2.0 * delta
        _, gminus, _ = calc.singlepoint(numbers, shifted)

        fd = ((gplus - gminus) / (2.0 * delta)).reshape(-1)
        assert np.abs(fd - hess[:, 0]).max() < 1.0e-5

    def test_rejects_periodic(self, sio2):
        numbers, positions, lattice, npbc = sio2
        calc = GFNFFCalculator(numbers, positions, lattice=lattice, npbc=npbc)
        with pytest.raises(NotImplementedError, match="periodic"):
            calc.hessian(numbers, positions)

    def test_step_argument_is_accepted(self, caffeine):
        numbers, positions = caffeine
        calc = GFNFFCalculator(numbers, positions)
        default, _, _ = calc.hessian(numbers, positions)
        tighter, _, _ = calc.hessian(numbers, positions, step=1.0e-3)
        # a different finite-difference step moves the numbers a little but
        # must not change the answer materially
        assert np.abs(default - tighter).max() < 1.0e-4


class TestCharges:
    """EEQ partial charges, exposed through c_gfnff_calculator_charges."""

    def test_sum_equals_total_charge(self, caffeine):
        numbers, positions = caffeine
        with GFNFFCalculator(numbers, positions) as calc:
            calc.singlepoint(numbers, positions)
            q = calc.charges()
        assert q.shape == (len(numbers),)
        assert q.sum() == pytest.approx(0.0, abs=1e-8)

    def test_sum_equals_nonzero_total_charge(self, caffeine):
        numbers, positions = caffeine
        with GFNFFCalculator(numbers, positions, charge=1) as calc:
            calc.singlepoint(numbers, positions)
            q = calc.charges()
        assert q.sum() == pytest.approx(1.0, abs=1e-8)

    def test_signs_follow_electronegativity(self, caffeine):
        """Caffeine's N and O must come out negative, its H positive.

        A shape-only check would pass on a buffer that was never written;
        this one fails if the charges are zeros or garbage.
        """
        numbers, positions = caffeine
        with GFNFFCalculator(numbers, positions) as calc:
            calc.singlepoint(numbers, positions)
            q = calc.charges()
        assert (q[numbers == 7] < 0).all()
        assert (q[numbers == 8] < 0).all()
        assert (q[numbers == 1] > 0).all()

    def test_tracks_geometry(self, caffeine):
        """The charges belong to the last singlepoint, not to init."""
        numbers, positions = caffeine
        with GFNFFCalculator(numbers, positions) as calc:
            calc.singlepoint(numbers, positions)
            first = calc.charges().copy()
            stretched = positions.copy()
            stretched[0, 0] += 0.3
            calc.singlepoint(numbers, stretched)
            second = calc.charges()
        assert np.abs(first - second).max() > 1e-6

    def test_requires_singlepoint(self, caffeine):
        numbers, positions = caffeine
        with GFNFFCalculator(numbers, positions) as calc:
            with pytest.raises(RuntimeError, match="singlepoint"):
                calc.charges()

    def test_pbc_charges(self, sio2):
        numbers, positions, lattice, npbc = sio2
        with GFNFFCalculator(
            numbers, positions, lattice=lattice, npbc=npbc
        ) as calc:
            calc.singlepoint(numbers, positions, lattice=lattice)
            q = calc.charges()
        assert q.sum() == pytest.approx(0.0, abs=1e-8)
        # Si donates to O in quartz
        assert (q[numbers == 14] > 0).all()
        assert (q[numbers == 8] < 0).all()


class TestVersion:
    """Parametrisation version selection."""

    def test_enum_matches_library(self):
        # the enum is built from the library, so this pins the values that the
        # C header documents
        assert int(Version.angewChem2020) == 1
        assert int(Version.angewChem2020_1) == 2
        assert int(Version.angewChem2020_2) == 3
        assert int(Version.mcgfnff2023) == 4
        assert int(Version.conformer2020) == 5
        assert int(Version.harmonic2020) == -1

    def test_accepts_name_int_and_enum(self, caffeine):
        numbers, positions = caffeine
        energies = []
        for spec in ("harmonic2020", -1, Version.harmonic2020):
            with GFNFFCalculator(numbers, positions, version=spec) as calc:
                energies.append(calc.singlepoint(numbers, positions)[0])
        assert energies[0] == energies[1] == energies[2]

    def test_unknown_name_rejected(self, caffeine):
        numbers, positions = caffeine
        with pytest.raises(ValueError, match="unknown GFN-FF version"):
            GFNFFCalculator(numbers, positions, version="no-such-version")

    def test_harmonic_changes_energy(self, caffeine):
        """harmonic2020 replaces the bond potential, so it must not agree."""
        numbers, positions = caffeine
        with GFNFFCalculator(numbers, positions) as calc:
            default = calc.singlepoint(numbers, positions)[0]
        with GFNFFCalculator(numbers, positions, version="harmonic2020") as calc:
            harmonic = calc.singlepoint(numbers, positions)[0]
        assert abs(harmonic - default) > 1e-6

    def test_conformer_matches_default_at_equilibrium(self, caffeine):
        """conformer2020 only replaces the dissociative tail of the bond term.

        Near equilibrium every bond is inside the region the two share, so the
        energy must not move.  Single threaded the two agree bit for bit; the
        tolerance is there because the OpenMP reductions sum in thread
        completion order, which costs a couple of ULP.
        """
        numbers, positions = caffeine
        with GFNFFCalculator(numbers, positions) as calc:
            default = calc.singlepoint(numbers, positions)[0]
        with GFNFFCalculator(numbers, positions, version="conformer2020") as calc:
            conformer = calc.singlepoint(numbers, positions)[0]
        assert abs(conformer - default) < 1e-12

    def test_conformer_bond_cannot_dissociate(self, caffeine):
        """Pull a bond well past the inflection point of the Gaussian well.

        The calculator is built once at the reference geometry, so the bond
        list is fixed and both versions are asked about the same bond.  The
        published well flattens out; the conformer variant keeps climbing.
        """
        numbers, positions = caffeine
        stretched = positions.copy()
        direction = positions[14] - positions[0]
        direction /= np.linalg.norm(direction)

        def curve(version):
            with GFNFFCalculator(numbers, positions, version=version) as calc:
                out = []
                for k in range(4):
                    stretched[14] = positions[14] + direction * 2.0 * k
                    out.append(calc.singlepoint(numbers, stretched)[0])
                return out

        ref = curve(None)
        con = curve("conformer2020")
        # the published term has run out of restoring force
        assert (ref[3] - ref[2]) / (ref[2] - ref[1]) < 0.1
        # the continuation is linear: equal steps cost equal energy
        assert 0.9 < (con[3] - con[2]) / (con[2] - con[1]) < 1.1
        assert con[3] - ref[3] > 0.1

    def test_harmonic_has_no_charges(self, caffeine):
        """It returns before the EEQ solve; the zeros must not be reported."""
        numbers, positions = caffeine
        with GFNFFCalculator(numbers, positions, version="harmonic2020") as calc:
            calc.singlepoint(numbers, positions)
            with pytest.raises(RuntimeError, match="no charges"):
                calc.charges()

    def test_mcgfnff_rejected_for_molecule(self, caffeine):
        numbers, positions = caffeine
        with pytest.raises(RuntimeError, match="molecular crystals"):
            GFNFFCalculator(numbers, positions, version="mcgfnff2023")

    def test_mcgfnff_accepted_for_crystal(self, sio2):
        numbers, positions, lattice, npbc = sio2
        with GFNFFCalculator(
            numbers, positions, lattice=lattice, npbc=npbc, version="mcgfnff2023"
        ) as calc:
            energy, _, _ = calc.singlepoint(numbers, positions, lattice=lattice)
        assert np.isfinite(energy)
        # mc adjusts the dispersion between fragments, so it must differ from
        # the default parametrisation
        with GFNFFCalculator(
            numbers, positions, lattice=lattice, npbc=npbc
        ) as calc:
            default, _, _ = calc.singlepoint(numbers, positions, lattice=lattice)
        assert abs(energy - default) > 1e-8


class TestParametrisation:
    """Parameter files supplied through the ``parametrisation`` argument."""

    def test_missing_file_rejected(self, caffeine):
        numbers, positions = caffeine
        with pytest.raises(FileNotFoundError):
            GFNFFCalculator(numbers, positions, parametrisation="/nope/x.toml")

    @pytest.mark.skipif(not toml_available(), reason="built without toml-f")
    def test_shipped_file_reproduces_internal(self, caffeine, shipped_param):
        """param/gfnff-angewChem2020_2.toml is a dump of the internal set."""
        numbers, positions = caffeine
        with GFNFFCalculator(numbers, positions) as calc:
            internal = calc.singlepoint(numbers, positions)[0]
        with GFNFFCalculator(
            numbers, positions, parametrisation=shipped_param
        ) as calc:
            from_file = calc.singlepoint(numbers, positions)[0]
        # not bit-exact: toml-f renders floats with 16 digits after the point
        # rather than 16 significant digits
        assert from_file == pytest.approx(internal, abs=1e-7)

    @pytest.mark.skipif(not toml_available(), reason="built without toml-f")
    def test_partial_overlay_applies(self, caffeine, tmp_path):
        """A file naming one key overlays the internal set for the rest."""
        numbers, positions = caffeine
        overlay = tmp_path / "tweak.toml"
        overlay.write_text("[generator]\nfringbo = 0.5\n")

        with GFNFFCalculator(numbers, positions) as calc:
            internal = calc.singlepoint(numbers, positions)[0]
        with GFNFFCalculator(
            numbers, positions, parametrisation=str(overlay)
        ) as calc:
            tweaked = calc.singlepoint(numbers, positions)[0]
        assert abs(tweaked - internal) > 1e-3


class TestBondMatrix:
    """Host-supplied molecular graph."""

    @staticmethod
    def _graph(positions, cutoff=3.2):
        """Adjacency from distances; on caffeine this is exactly the bonds."""
        n = len(positions)
        d = np.linalg.norm(positions[:, None] - positions[None, :], axis=-1)
        bm = ((d < cutoff) & (d > 0.0)).astype(np.int32)
        return bm

    def test_graph_reproduces_perception(self, caffeine):
        """The graph GFN-FF would find, handed back, must change nothing."""
        numbers, positions = caffeine
        with GFNFFCalculator(numbers, positions) as calc:
            default = calc.singlepoint(numbers, positions)[0]
        bm = self._graph(positions)
        with GFNFFCalculator(numbers, positions, bond_matrix=bm) as calc:
            supplied = calc.singlepoint(numbers, positions)[0]
        assert abs(supplied - default) < 1e-12

    def test_bond_orders_are_not_consumed(self, caffeine):
        """Only the zero/nonzero pattern has an effect, as documented."""
        numbers, positions = caffeine
        bm = self._graph(positions)
        with GFNFFCalculator(numbers, positions, bond_matrix=bm) as calc:
            single = calc.singlepoint(numbers, positions)[0]
        with GFNFFCalculator(numbers, positions, bond_matrix=bm * 3) as calc:
            tripled = calc.singlepoint(numbers, positions)[0]
        assert abs(tripled - single) < 1e-12

    def test_wrong_shape_rejected(self, caffeine):
        numbers, positions = caffeine
        with pytest.raises(ValueError, match="bond_matrix must have shape"):
            GFNFFCalculator(numbers, positions, bond_matrix=np.zeros((3, 3), int))

    def test_asymmetric_rejected(self, caffeine):
        numbers, positions = caffeine
        bm = self._graph(positions)
        bm[0, 1] = 1
        bm[1, 0] = 0
        with pytest.raises(RuntimeError):
            GFNFFCalculator(numbers, positions, bond_matrix=bm)

    def test_self_bond_rejected(self, caffeine):
        numbers, positions = caffeine
        bm = self._graph(positions)
        bm[0, 0] = 1
        with pytest.raises(RuntimeError):
            GFNFFCalculator(numbers, positions, bond_matrix=bm)
