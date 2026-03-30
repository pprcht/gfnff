<div align="center">

<h1>GFN-FF</h1>
<h3>A general force field for elements <i>Z</i> = 1–103</h3>

![build status](https://github.com/pprcht/gfnff/actions/workflows/build-and-test.yml/badge.svg)
[![License: LGPL v3](https://img.shields.io/badge/license-LGPL_v3-coral.svg)](https://www.gnu.org/licenses/lgpl-3.0)

</div>

This repository provides a standalone library implementation of the **GFN-FF** method by **S.Spicher** and **S.Grimme**,
adapted from the [`xtb`](https://github.com/grimme-lab/xtb) code (most recently at commit [`6d44803`](https://github.com/grimme-lab/xtb/commit/6d44803)
and validated against that version's results).
The primary purpose is to serve as a linkable dependency for other Fortran, C, and C++ projects.
From this point forward, development may diverge from the upstream `xtb` implementation.

---

## Method

GFN-FF (*Geometries, Frequencies, Non-covalent interactions Force-Field*) is a
completely automated, topology-based force field for fast structure optimisations
and non-covalent interaction energies across the periodic table (Z = 1–103).
The topology and parametrisation are derived entirely from the input geometry,
without user-defined atom types or connectivity.

The following features are available and documented in the associated publications:

- **Molecular GFN-FF** — a generic, partially polarisable force field covering organic,
  organometallic, and biochemical systems
  (S. Spicher, S. Grimme, *Angew. Chem. Int. Ed.* **2020**, 59, 15665.
  [doi:10.1002/anie.202004239](https://doi.org/10.1002/anie.202004239))

- **Periodic boundary conditions / molecular crystals** — adjusted non-covalent
  interactions for lattice energy predictions and unit-cell optimisations of molecular crystals
  (S. Grimme, T. Rose, *Z. Naturforsch. B* **2024**, 79, 191.
  [doi:10.1515/znb-2023-0088](https://doi.org/10.1515/znb-2023-0088))

- **Lanthanide and actinide extension** — reparametrised f-element treatment enabling
  MD simulations and geometry optimisations for large lanthanide- and actinide-containing
  systems
  (T. Rose, M. Bursch, J.-M. Mewes, S. Grimme, *Inorg. Chem.* **2024**.
  [doi:10.1021/acs.inorgchem.4c03215](https://doi.org/10.1021/acs.inorgchem.4c03215))

---

## Building from source

The library requires a Fortran and C compiler (e.g. `gfortran`/`gcc`),
LAPACK/BLAS (e.g. OpenBLAS), and optionally OpenMP.
Both CMake (≥ 3.21) and Meson (≥ 0.59) build systems are supported.

<table>
<tr>
<th>CMake</th>
<th>Meson</th>
</tr>
<tr>
<td>

```bash
cmake -B _build
cmake --build _build
```

To run the test suite:

```bash
cmake -B _build -DWITH_TESTS=ON
cmake --build _build
ctest --test-dir _build
```

</td>
<td>

```bash
meson setup _build
ninja -C _build
```

To run the test suite:

```bash
meson setup _build -Dtests=true
ninja -C _build test
```

</td>
</tr>
</table>

The compiled library (`libgfnff.a` by default) is placed in the build directory
and can be linked into any downstream project.

---

## Library usage

The interface is exposed through the `gfnff_interface` Fortran module and the
`gfnff_interface_c.h` C header (located in `include/`).
Two steps are required: initialise the calculator (topology setup) and call the
singlepoint routine. The initialisation is typically the more expensive step;
once complete, singlepoint evaluations can be called repeatedly on the same
calculator object.

<table>
<tr><td>
<details>
<summary><b>Fortran</b></summary>

```fortran
use iso_fortran_env, only: real64
use gfnff_interface

type(gfnff_data) :: calc
integer  :: nat, ichrg, io
integer,  allocatable :: at(:)
real(real64), allocatable :: xyz(:,:), gradient(:,:)
real(real64) :: energy, sigma(3,3)

! ... populate nat, at, xyz, ichrg ...

call calc%init(nat, at, xyz, ichrg=ichrg, iostat=io)

call calc%singlepoint(nat, at, xyz, energy, gradient, iostat=io, sigma=sigma)

call calc%deallocate()
```

All coordinates are in Bohr; the energy is in Hartree, the gradient in Eh/Bohr,
and `sigma` (3×3) is the stress tensor in Hartree (zero for non-periodic systems).
Full working example: [`app/main.F90`](app/main.F90).

</details>
</td></tr>
<tr><td>
<details>
<summary><b>C</b></summary>

```c
#include "gfnff_interface_c.h"

double sigma[3][3];   /* stress tensor (Hartree); zero for non-PBC */

c_gfnff_calculator calc =
    c_gfnff_calculator_init(nat, at, xyz, ichrg, printlevel, solvent);

c_gfnff_calculator_singlepoint(&calc, nat, at, xyz, &energy, gradient,
                               sigma, &iostat);

c_gfnff_calculator_deallocate(&calc);
```

Full working example: [`test/main.c`](test/main.c).

</details>
</td></tr>
<tr><td>
<details>
<summary><b>C++</b></summary>

```cpp
#include "gfnff_interface_c.h"

double sigma[3][3];   // stress tensor (Hartree); zero for non-PBC

c_gfnff_calculator calc =
    c_gfnff_calculator_init(nat, at, xyz, ichrg, printlevel, solvent);

c_gfnff_calculator_singlepoint(&calc, nat, at, xyz, &energy, gradient,
                               sigma, &iostat);

c_gfnff_calculator_deallocate(&calc);
```

Full working example: [`test/main.cpp`](test/main.cpp).

</details>
</td></tr>
<tr><td>
<details>
<summary><b>Integrating as a CMake subproject</b></summary>

Add the repository as a subdirectory and link against the exported target:

```cmake
add_subdirectory(gfnff)
target_link_libraries(my_target PRIVATE gfnff)
```

</details>
</td></tr>
<tr><td>
<details>
<summary><b>Integrating as a Meson subproject</b></summary>

Place the repository under `subprojects/gfnff/` and wrap it:

```meson
gfnff_dep = dependency('gfnff', fallback: ['gfnff', 'gfnff_dep'])
```

</details>
</td></tr>
</table>

---

## Periodic boundary conditions

PBC support is available via `c_gfnff_calculator_init_pbc` on the C/C++ side
and via an optional `lattice` argument to `calc%init` in Fortran.
See the PBC sections in [`test/main.c`](test/main.c) and [`test/main.cpp`](test/main.cpp)
for worked examples.

---

## Python bindings

Python bindings are provided through a `ctypes`-based interface.
The shared library is bundled into a binary wheel, so no Fortran or C compiler
is needed at install time.

### Installation

Binary wheels for Linux (x86\_64) and macOS (x86\_64 / arm64) are published on PyPI:

```bash
pip install gfnff          # library only
pip install "gfnff[ase]"   # + ASE (enables the CLI and the ASE calculator)
```

### Building from source

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

### Command-line interface

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

---

### Low-level API (`GFNFFCalculator`)

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

### ASE Calculator (`GFNFF`)

`GFNFF` is a fully compatible [ASE `Calculator`](https://wiki.fysik.dtu.dk/ase/ase/calculators/calculators.html).
It handles unit conversion automatically (Å ↔ Bohr, eV ↔ Hartree).
Implemented properties: **energy**, **forces**, **stress**.

```python
from ase.build import molecule
from gfnff import GFNFF

atoms = molecule("caffeine")
atoms.calc = GFNFF()

energy = atoms.get_potential_energy()   # eV
forces = atoms.get_forces()             # eV / Å, shape (nat, 3)
stress = atoms.get_stress()             # eV / Å³, Voigt [xx,yy,zz,yz,xz,xy]; zero for non-PBC
```

Periodic systems work the same way — provide an `atoms` object with `cell` and `pbc` set:

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

Additional options:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `charge` | `0` | Total charge. Also reads `atoms.info["charge"]` (takes precedence). |
| `solvent` | `""` | Implicit solvent name: `"h2o"`, `"acetone"`, `"chcl3"`, … (molecular systems only) |
| `printlevel` | `0` | Fortran output verbosity (0 = silent, 3 = verbose). |

### Running the tests

```bash
pip install "gfnff[test]"
pytest python/tests/
```

---

## License

This project is licensed (as the original `xtb` code) under the **GNU Lesser General Public License v3** or later.
See [`LICENSE`](LICENSE) for details.
