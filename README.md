
# GFN-FF

This repository contains a standalone implementation of
the GFN-FF method by S.Spicher and S.Grimme (<https://doi.org/10.1002/anie.202004239>), adapted from the [`xtb` code](https://github.com/grimme-lab/xtb).

The default CMake and meson builds compile a statically linked library (`libgfnff.a`) that can be linked in other projects.

[`main.f90`](testprog/main.f90) in [`testprog/`](testprog/) demonstrates the in-code usage.


---

## Building the Project

Make sure you have the following dependencies installed:

- CMake and `make`, or meson and ninja build systems
- Fortran and C compilers (e.g., `gfortran`/`gcc` or `ifort`/`icc`)

### Instructions

Follow these steps to build the project:

1. Create a build directory and navigate to it
   ```bash
   mkdir _build
   cd _build
   ```

2. Export the compilers (here for example `ifort`/`icc`) and depending on your chosen build system   set up the build:
   - generate the build files using CMake:
     ```bash
     FC=ifort CC=icc cmake ..
     ```
   - generate the build files using meson:
     ```bash
     FC=ifort CC=icc meson ..
     ```
   I you wish to build the test-binary, add `-Dbuild_exe=true` to either the `cmake` or `meson`      setup command.


3. Depending on your chosen build system, build the project. If you have multiple cores/processors,  you can speed up the build process by specifying the number of cores to use with the `-j` option.    For example, to use 4 cores:
   - With CMake/`make`:
     ```shell
     make -j4
     ```
   - With meson/`ninja`:
     ```shell
     ninja -j4
     ```
### Cleaning the Build

To clean the build files, simply delete the `build` directory:

```shell
rm -rf _build
```


