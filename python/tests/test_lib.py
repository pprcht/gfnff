# This file is part of gfnff.
# SPDX-Identifier: LGPL-3.0-or-later
"""
Unit tests for the low-level GFNFFCalculator (Bohr / Hartree units).
These mirror the test cases in test/main.c.
"""

import numpy as np
import pytest

from gfnff import GFNFFCalculator

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
