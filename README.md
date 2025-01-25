
# GFN-FF - A general force-field for elements *Z*=1-86

This repository contains a standalone implementation of
the **GFN-FF** method by **S.Spicher** and **S.Grimme** (<https://doi.org/10.1002/anie.202004239>), and was adapted from the [`xtb` code](https://github.com/grimme-lab/xtb).

The default CMake build compiles a statically linked library (`libgfnff.a`) that can be linked in other projects.

[`main.f90`](app/main.f90) in [`app/`](app/) demonstrates the (Fortran) in-code usage.


---

### Instructions (building from source)

Make sure you have the following dependencies installed:

- CMake and `make`
- Fortran and C compilers (e.g., `gfortran`/`gcc`)

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

`libgfnff.a` (and a `gfnff` app, when `-Dbuild_exe=true`) can now be found in this (`_build`) directory.

---

