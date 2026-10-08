# Library usage

The interface is exposed through the `gfnff_interface` Fortran module and the
`gfnff_interface_c.h` C header (located in `include/`).
Two steps are required: initialise the calculator (topology setup) and call the
singlepoint routine. The initialisation is typically the more expensive step;
once complete, singlepoint evaluations can be called repeatedly on the same
calculator object.

All coordinates are in Bohr; the energy is in Hartree, the gradient in Eh/Bohr,
and `sigma` (3×3) is the stress tensor in Hartree (zero for non-periodic systems).

## Fortran

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

Full working example: [`app/main.F90`](../app/main.F90).

## C and C++

The same header serves both languages.

```c
#include "gfnff_interface_c.h"

double sigma[3][3];   /* stress tensor (Hartree); zero for non-PBC */

c_gfnff_calculator calc =
    c_gfnff_calculator_init(nat, at, xyz, ichrg, printlevel, solvent);

c_gfnff_calculator_singlepoint(&calc, nat, at, xyz, &energy, gradient,
                               sigma, &iostat);

c_gfnff_calculator_deallocate(&calc);
```

Full working examples: [`test/main.c`](../test/main.c) and
[`test/main.cpp`](../test/main.cpp).

## Periodic boundary conditions

PBC support is available via `c_gfnff_calculator_init_pbc` on the C/C++ side
and via an optional `lattice` argument to `calc%init` in Fortran.
See the PBC sections in [`test/main.c`](../test/main.c) and
[`test/main.cpp`](../test/main.cpp) for worked examples.

## Integrating into another project

**CMake.** Add the repository as a subdirectory and link against the exported target:

```cmake
add_subdirectory(gfnff)
target_link_libraries(my_target PRIVATE gfnff)
```

**Meson.** Place the repository under `subprojects/gfnff/` and wrap it:

```meson
gfnff_dep = dependency('gfnff', fallback: ['gfnff', 'gfnff_dep'])
```
