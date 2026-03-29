# Python Bindings — Design Notes

This document explains the key decisions made when adding Python bindings to
GFN-FF.  The goal is to make the rationale reproducible so that future
maintainers can judge edge cases without having to rediscover the reasoning.

---

## ctypes instead of pybind11

The C API (`include/gfnff_interface_c.h`) already exposes a clean, stable,
5-function interface.  Two options were considered:

| | ctypes | pybind11 |
|---|---|---|
| Extra compiled code | none — loads the shared library at import time | a C++ extension module must be compiled |
| User toolchain | just Python + pip | Fortran + C++ compiler at install time |
| API duplication | zero — the C API is used as-is | another wrapper layer around the existing C wrapper |
| Type safety | runtime checks via `ndpointer` | compile-time via C++ templates |
| Maintenance | one layer (Python ↔ C API) | two layers (Python ↔ C++ ↔ C API) |

**Decision**: ctypes.

pybind11's main strength — wrapping C++ classes with compile-time type safety —
is wasted here because the Fortran objects are already hidden behind an opaque
`void *` pointer.  The only real gain would be a slightly more Pythonic API, but
that is achieved equally well in pure Python on top of ctypes.

The commented-out pybind11 block that existed in `CMakeLists.txt` was replaced
with the shared-library target described below.

---

## Static library + separate shared library

The existing build produced a static archive (`libgfnff.a`) with
`POSITION_INDEPENDENT_CODE=TRUE` — ready for downstream Fortran and CMake
consumers.  ctypes requires a shared library, so a parallel `gfnff-shared`
target was added, guarded by `-DPYTHON_BINDINGS=ON`.

**Why not just switch to a shared library globally?**  The static library is the
right default for downstream CMake/Meson subproject use: it avoids runtime path
issues and produces self-contained executables.  The shared library is only
needed for the Python wheel.

**Why a separate Fortran module directory?**  Both targets compile the same
source files.  If they shared `${CMAKE_CURRENT_BINARY_DIR}/include` for their
`.mod` files, parallel builds could race.  The shared target uses
`include-shared/` instead.

---

## scikit-build-core as the wheel build backend

`scikit-build-core` (the modern replacement for the deprecated `scikit-build`)
drives CMake from inside the PEP 517 build process.  The relevant flow is:

1. `pip install gfnff` (or `python -m build`) invokes scikit-build-core.
2. scikit-build-core calls `cmake -DPYTHON_BINDINGS=ON …` and `cmake --build`.
3. The CMake `POST_BUILD` command copies `libgfnff.so` into `python/gfnff/`
   (see next section).
4. scikit-build-core packages `python/gfnff/` via `wheel.packages`, bundling
   the library alongside the Python sources.
5. At import time `_lib.py` finds `libgfnff.so*` next to `__init__.py` and
   loads it with `ctypes.CDLL`.

Plain `setuptools` was not used because it has no native support for CMake
builds; it would require a custom `build_ext` subclass or a subprocess call,
which is more fragile.

---

## POST_BUILD copy instead of CMake install() for the shared library

**Problem**: ctypes needs the shared library to live next to the Python package
files so `_lib.py` can find it with a simple `glob("libgfnff.so*")`.  The
natural CMake solution is an `install(TARGETS … DESTINATION gfnff)` rule, which
places the library in the right spot during `cmake --install`.  This works for
wheel builds but breaks the development workflow: editable installs (`pip
install -e .`) never run `cmake --install`, so the `.so` stays in the build
directory and every Python invocation requires `LD_LIBRARY_PATH=_build_py/`.

**Solution**: replace the `install()` rule with a `POST_BUILD` custom command:

```cmake
add_custom_command(
  TARGET "${PROJECT_NAME}-shared" POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy_if_different
    $<TARGET_FILE:${PROJECT_NAME}-shared>
    "${CMAKE_CURRENT_SOURCE_DIR}/python/gfnff/$<TARGET_FILE_NAME:${PROJECT_NAME}-shared>"
)
```

`POST_BUILD` runs after every successful `cmake --build`, regardless of whether
an install step follows.  The library lands in `python/gfnff/` immediately,
where `_lib.py` finds it without any environment variable.

**Why keep `install()` as well?**  scikit-build-core's `wheel.packages` collects
Python source files; it does not reliably include compiled binaries like `.so`.
The CMake `install(TARGETS gfnff-shared LIBRARY DESTINATION gfnff)` step is
what actually gets the shared library into the wheel.  The `POST_BUILD` copy is
only for local editable installs where `cmake --install` is never run.

**`python/gfnff/libgfnff*` is gitignored** because it is a build artefact
generated from sources already tracked in the repository.

---

## Array memory layout (row-major / column-major)

NumPy's default layout is C row-major.  A `(nat, 3)` float64 array stores the
three coordinates of atom 0 contiguously, then atom 1, etc. — exactly the same
as the C declaration `double xyz[nat][3]`.

The Fortran side stores arrays column-major: `xyz(3, nat)` with coordinates of
atom 0 in column 1.  The `C_interface.f90` wrapper already handles this
transpose internally via `c_f_pointer` with the reversed shape.

**Consequence**: no explicit transpose is needed in the Python code.  Passing a
C-contiguous `(nat, 3)` NumPy array directly to `c_gfnff_calculator_singlepoint`
is correct.  The `ndpointer(flags="C_CONTIGUOUS")` validator in `_lib.py`
enforces this; `np.ascontiguousarray()` is called before every ctypes call to
guarantee it even if the caller passed a sliced or Fortran-order array.

---

## Topology re-initialisation heuristic (ASE calculator)

GFN-FF separates topology setup (expensive: bond detection, ring detection,
parameter assignment) from the energy/gradient evaluation (cheap: force
evaluation on fixed topology).  The Fortran data structure holds the topology
between singlepoint calls.

The ASE `Calculator.calculate` method receives a `system_changes` list that
describes what changed since the last call.  The heuristic used:

- **Re-initialise** when `numbers` changed (different element count or types) or
  when `cell`/`pbc` appears in `system_changes` (periodic boundary conditions
  changed).
- **Reuse topology** for pure geometry updates (`positions` in
  `system_changes` only).

This matches the intended usage of GFN-FF: the topology is geometry-dependent
only through the initial bond-detection step, which uses distance thresholds
rather than exact positions.  Small displacements during an optimisation do not
invalidate the topology.

---

## Stress tensor (deferred)

The Fortran `gfnff_singlepoint` subroutine already computes the stress tensor
(`sigma(3,3)`, in Eh).  However, `c_gfnff_calculator_singlepoint` in
`C_interface.f90` does not expose it as an output argument.

Extending the C API requires:
1. Adding `double sigma[3][3]` (or a NULL-pointer variant) to the C signature.
2. Passing it through the Fortran wrapper.
3. Converting to Voigt notation in the ASE calculator using
   `ase.stress.full_3x3_to_voigt_6_stress` and dividing by cell volume.

Until then the ASE calculator raises `PropertyNotImplementedError` when
`"stress"` is requested, which is the correct ASE pattern for missing
properties.

---

## manylinux wheel static linking

Linux wheels must pass `auditwheel`, which rejects wheels that dynamically link
against libraries outside the manylinux allowlist.  The Fortran runtime
(`libgfortran`, `libquadmath`, `libgomp`) and LAPACK (OpenBLAS) are not on that
list.

The CI workflow passes:
```
-DCMAKE_Fortran_FLAGS=-static-libgfortran
-DBLA_STATIC=ON
```
to CMake so these dependencies are baked into `libgfnff.so` itself.
`auditwheel repair` (run automatically by cibuildwheel) then bundles any
remaining allowed shared libraries (e.g. `libc`, `libpthread`) into the wheel.

---

## Conda-forge

Conda-forge packages are maintained as independent feedstock repositories
(`https://github.com/conda-forge/<name>-feedstock`).  The recipe lives there,
not in this repository.  A feedstock should be opened once the first stable
version is published on PyPI.
