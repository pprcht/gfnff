# Python bindings

Python bindings are provided through a `ctypes`-based interface.
The shared library is bundled into a binary wheel, so no Fortran or C compiler
is needed at install time.

## Installation

Binary wheels for Linux (x86\_64, aarch64) and macOS (arm64) are published on PyPI:

```bash
pip install gfnff          # library only
pip install "gfnff[ase]"   # + ASE (enables the CLI and the ASE calculator)
```

## Building from source

The source build compiles the Fortran library on your machine.
The following **system packages** must be present before running pip:

| Dependency | Example (Debian/Ubuntu) | Example (Fedora/RHEL) | Example (macOS) |
|---|---|---|---|
| Fortran compiler | `apt install gfortran` | `dnf install gcc-gfortran` | `brew install gcc` |
| LAPACK + BLAS | `apt install libopenblas-dev` | `dnf install openblas-devel` | `brew install openblas` |
| CMake ≥ 3.21 | installed by pip automatically | ← | ← |

Once those are in place:

```bash
pip install ".[ase]"           # from a checkout
pip install "gfnff[ase]" --no-binary gfnff   # force source build from PyPI
```

## Command-line interface

Installing `gfnff[ase]` places a `gfnff` executable on your PATH.

```
gfnff <input> [options]
```

The input file is read by ASE, so any format it supports works (xyz, extxyz, POSCAR, cif, …).

**Singlepoint** (default):

```bash
gfnff molecule.xyz
gfnff molecule.xyz --chrg -1
gfnff molecule.xyz --alpb h2o        # implicit solvation (--solv is an alias)
```

**Geometry optimisation** (L-BFGS via ASE, cell fixed, writes `gfnff.log.extxyz`):

```bash
gfnff molecule.xyz --opt
gfnff molecule.xyz --opt --fmax 0.05          # looser convergence, eV/Å
gfnff molecule.xyz --opt --outfile path.xyz   # custom trajectory file
gfnff molecule.xyz --opt --alpb acetone       # optimise in solvent
```

**Variable-cell optimisation** (L-BFGS + `ExpCellFilter`, periodic systems only):

```bash
gfnff crystal.cif --optcell
gfnff crystal.cif --optcell --fmax 0.01
```

Full option list: `gfnff --help`

The trajectory file (`gfnff.log.extxyz`) stores energy and forces in each
frame header, compatible with ASE's `ase gui`.

## Low-level API (`GFNFFCalculator`)

`GFNFFCalculator` mirrors the C API one-to-one.
All quantities use the same units as the library itself: **Bohr** for coordinates and lattice, **Hartree** for energy, **Eh/Bohr** for gradients, and **Hartree** for the stress tensor.

```python
import numpy as np
from gfnff import GFNFFCalculator

# Atomic numbers and coordinates in Bohr
numbers = np.array([6, 8, 1, 1], dtype=np.int32)   # CO + 2 H
positions = np.array([[0, 0, 0], [2.1, 0, 0],
                      [-1.0, 0, 0], [3.1, 0, 0]], dtype=np.float64)

with GFNFFCalculator(numbers, positions, charge=0, printlevel=0) as calc:
    energy, gradient, sigma = calc.singlepoint(numbers, positions)
    print(f"Energy: {energy:.6f} Eh")
    print(f"Gradient shape: {gradient.shape}")  # (nat, 3)
    print(f"Stress tensor:\n{sigma}")            # (3, 3), Hartree; zero for non-PBC
```

Periodic systems use a separate initialiser:

```python
calc = GFNFFCalculator(
    numbers, positions,
    lattice=lattice_bohr,   # shape (3, 3), rows are lattice vectors
    npbc=3,
)
```

Partial charges are those of the last singlepoint, where they are obtained as
part of the energy evaluation:

```python
with GFNFFCalculator(numbers, positions) as calc:
    calc.singlepoint(numbers, positions)
    q = calc.charges()      # (nat,), in e; sums to the total charge
```

Force-field versions, custom parameter files and user-supplied molecular graphs
are described in [parametrisation.md](parametrisation.md).

## ASE Calculator (`GFNFF`)

`GFNFF` is a fully compatible [ASE `Calculator`](https://wiki.fysik.dtu.dk/ase/ase/calculators/calculators.html).
It handles unit conversion automatically (Å ↔ Bohr, eV ↔ Hartree).
Implemented properties: **energy**, **forces**, **stress**, **charges**.

```python
from ase.build import molecule
from gfnff import GFNFF

atoms = molecule("caffeine")
atoms.calc = GFNFF()

energy = atoms.get_potential_energy()   # eV
forces = atoms.get_forces()             # eV / Å, shape (nat, 3)
stress = atoms.get_stress()             # eV / Å³, Voigt [xx,yy,zz,yz,xz,xy]; zero for non-PBC
charges = atoms.get_charges()           # e, shape (nat,)
```

Periodic systems work the same way; provide an `atoms` object with `cell` and `pbc` set:

```python
from ase.io import read
from gfnff import GFNFF

atoms = read("quartz.cif")
atoms.calc = GFNFF()
print(atoms.get_potential_energy())  # eV / unit cell
print(atoms.get_stress())            # eV / Å³, Voigt
```

Variable-cell relaxation via ASE's `ExpCellFilter`:

```python
from ase.filters import ExpCellFilter
from ase.optimize import LBFGS

opt = LBFGS(ExpCellFilter(atoms))
opt.run(fmax=0.01)
```

### Hessian and vibrations

The Hessian is an explicit method rather than one of `implemented_properties`,
so it is never computed as a side effect of an MD or optimisation step:

```python
hessian = atoms.calc.get_hessian(atoms)   # (3N, 3N), eV / Å²

vib = atoms.calc.get_vibrations(atoms)    # ase.vibrations.VibrationsData
vib.get_energies()                        # frequencies, normal modes,
vib.get_modes()                           # thermochemistry, all from ASE
```

`get_vibrations()` passes the analytic Hessian directly to ASE, which avoids the 6N
displaced singlepoints that `ase.vibrations.Vibrations` would otherwise run. Results
for an unchanged geometry are cached, so asking for frequencies and modes does
not evaluate it twice. Non-periodic systems only.

Additional options:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `charge` | `0` | Total charge. Also reads `atoms.info["charge"]` (takes precedence). |
| `solvent` | `""` | Implicit solvent name: `"h2o"`, `"acetone"`, `"chcl3"`, … (molecular systems only) |
| `printlevel` | `0` | Fortran output verbosity (0 = silent, 3 = verbose). |
| `version` | `None` | Force-field version; see [parametrisation.md](parametrisation.md). |
| `parametrisation` | `None` | Path to a parameter file, overlaid on that version. |

Changing any of these with `calc.set(...)` rebuilds the underlying force field
rather than reusing the cached result.

## Running the tests

```bash
pip install "gfnff[test]"
pytest python/tests/
```
