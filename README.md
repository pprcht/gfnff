
# GFN-FF - A general force-field for elements *Z*=1-86

<div align="center">

![build status](https://github.com/pprcht/gfnff/actions/workflows/build-and-test.yml/badge.svg)
[![License: LGPL v3](https://img.shields.io/badge/license-LGPL_v3-coral.svg)](https://www.gnu.org/licenses/lgpl-3.0)

</div>


This repository contains a standalone implementation of
the **GFN-FF** method by **S.Spicher** and **S.Grimme** (<https://doi.org/10.1002/anie.202004239>), and was adapted from the [`xtb` code](https://github.com/grimme-lab/xtb).

The default CMake build compiles a statically linked library (`libgfnff.a`) that can be linked in other Fortran, C and C++ projects.




---

### Instructions (building from source)

Make sure you have the following dependencies installed:

- CMake and `make`
- Fortran and C compilers (e.g., `gfortran`/`gcc`), including LAPACK and OpenMP libraries (e.g. openblas)

Follow these steps to build the project assuming you are currently within the cloned repository:

1. Create a build directory and navigate to it
   ```bash
   mkdir _build
   cd _build
   ```

2. Run CMake to set up the build:
   - either directly:
     ```bash
     cmake ..
     ```
   - or if you wish to build the minimal example app:
     ```bash
     cmake .. -Dbuild_exe=true
     ```

3. Build the project and run the testsuite:
   ```
   make
   make test
   ``` 

`libgfnff.a` (and a `gfnff` app, when `-Dbuild_exe=true` was used in the CMake setup) can now be found in this (`_build`) directory.


---

### In-code instructions

[`main.f90`](app/main.f90) in [`app/`](app/) demonstrates the (Fortran) in-code usage.
[`main.c`](test/main.c) and [`main.cpp`](test/main.cpp) in [`test/`](test/) demonstrate the C and C++ in-code usage.

Two steps are required for using the GFN-FF library: 1. Initializing the calculator class (which stores topology and parametrization) and 2. calling the energy+gradient routine from it. A minimal Fortran example would look something like this:

```Fortran
use gfnff_interface
type(gfnff_data) :: calculator   !> model storage
integer :: nat                   !> number of atoms in the system 
integer :: at(nat)               !> atomic number for each atom
real(kind=real64) :: xyz(3,nat)  !> atomic coordinates in BOHR 
integer :: ichrg                 !> molecular charge
real(kind=real64) :: energy      !> energy in Hartree
real(kind=real64) :: grad(3,nat) !> gradient in Hartree/Bohr
integer :: io                    !> output status   

!> Read-in/define the system 
[...]

!> Initialize the calculator
call calculator%init(nat,at,xyz,ichrg=ichrg)

!> Calculate the energy and gradient routine
call calculator%singlepoint(nat,at,xyz,energy,gradient,iostatus=io)

!> Destroy calculator (once you are done)
call calculator%deallocate()

```

The calculator initialization is often more expensive than the actual energy evaluation due to the topology setup. However, once the calculator has been initialized, singlepoint energies can be called repeatedly.


