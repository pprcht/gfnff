<div align="center">

<h1>GFN-FF</h1>
<h3>A general force field for elements <i>Z</i> = 1–103</h3>

![build status](https://github.com/pprcht/gfnff/actions/workflows/build-and-test.yml/badge.svg)
[![License: LGPL v3](https://img.shields.io/badge/license-LGPL_v3-coral.svg)](https://www.gnu.org/licenses/lgpl-3.0)

</div>

This repository provides a standalone library implementation of the **GFN-FF** method,
adapted from the [`xtb`](https://github.com/grimme-lab/xtb) code (most recently at commit `6d44803`
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

To also build the standalone app:

```bash
cmake -B _build -Dbuild_exe=true
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

To also build the standalone app:

```bash
meson setup _build -Dbuild_exe=true
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

### Fortran

```fortran
use iso_fortran_env, only: real64
use gfnff_interface

type(gfnff_data) :: calc
integer  :: nat, ichrg, io
integer,  allocatable :: at(:)
real(real64), allocatable :: xyz(:,:), gradient(:,:)
real(real64) :: energy

! ... populate nat, at, xyz, ichrg ...

call calc%init(nat, at, xyz, ichrg=ichrg, iostat=io)

call calc%singlepoint(nat, at, xyz, energy, gradient, iostat=io)

call calc%deallocate()
```

All coordinates are in Bohr; the energy is in Hartree and the gradient in Eh/Bohr.

### C

```c
#include "gfnff_interface_c.h"

c_gfnff_calculator calc =
    c_gfnff_calculator_init(nat, at, xyz, ichrg, printlevel, solvent);

c_gfnff_calculator_singlepoint(&calc, nat, at, xyz, &energy, gradient, &iostat);

c_gfnff_calculator_deallocate(&calc);
```

### C++

```cpp
#include "gfnff_interface_c.h"

c_gfnff_calculator calc =
    c_gfnff_calculator_init(nat, at, xyz, ichrg, printlevel, solvent);

c_gfnff_calculator_singlepoint(&calc, nat, at, xyz, &energy, gradient, &iostat);

c_gfnff_calculator_deallocate(&calc);
```

Full working examples are available in [`app/main.F90`](app/main.F90) (Fortran),
[`test/main.c`](test/main.c) (C), and [`test/main.cpp`](test/main.cpp) (C++).

### Integrating as a CMake subproject

Add the repository as a subdirectory and link against the exported target:

```cmake
add_subdirectory(gfnff)
target_link_libraries(my_target PRIVATE gfnff)
```

### Integrating as a Meson subproject

Place the repository under `subprojects/gfnff/` and wrap it:

```meson
gfnff_dep = dependency('gfnff', fallback: ['gfnff', 'gfnff_dep'])
```

---

## Periodic boundary conditions

PBC support is available via `c_gfnff_calculator_init_pbc` on the C/C++ side
and via an optional `lattice` argument to `calc%init` in Fortran.
See the PBC sections in [`test/main.c`](test/main.c) and [`test/main.cpp`](test/main.cpp)
for worked examples.

---

## License

This project is licensed (as the original `xtb` code) under the **GNU Lesser General Public License v3** or later.
See [`LICENSE`](LICENSE) for details.
