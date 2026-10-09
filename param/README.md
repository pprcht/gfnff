# Parameter sets

TOML parametrisations. Point at one with `gfnff_data%parametrisation` (Fortran),
the `parametrisation` argument of `c_gfnff_calculator_init_ex` (C), or
`parametrisation=` on `GFNFFCalculator` / the ASE `GFNFF` calculator (Python).

A file carries both halves of a parametrisation:

- `[elements]` — the eleven per-element arrays for Z = 1..86 (`chi`, `gam`,
  `cnf`, `alp`, `bond`, `repa`, `repan`, `angl`, `angl2`, `tors`, `tors2`)
- `[generator]` — the global constants: ring prefactors, bond steepness,
  Hückel iteration, metal shifts, and the rest of `TGFFGenerator`

Physical constants (`en`, `rad`, `rcov`, `metal`, `group`, `normcn`, `repz`)
are deliberately absent. They are shared by every parametrisation and come
from `gfnff_param_tables`; round-tripping them would only invite drift.

## Partial files

A file need name only the keys it changes — the rest keep the values of the
internal parametrisation for the selected version. This is an overlay, so a
three-line file is a valid parameter set:

```toml
[generator]
fringbo = 0.5
```

## Writing one

```fortran
call gfnff_write_parametrisation(calc, 'my-set.toml', 'my label', iostat, errmsg)
```

`gfnff-angewChem2020_2.toml` is the library default, written by exactly that
call, and is the file to copy when starting a new set.

## Caveats

- Values round-trip to ~1e-16 relative, not bit-exactly: toml-f renders a
  float below 1e3 with sixteen digits *after the point* rather than sixteen
  significant digits. The effect on an energy is one unit in the last place.
- Parameter I/O must not be called from inside an OpenMP parallel region;
  see the note at the top of `src/core/param_io.F90`.
- Without toml-f the reader and writer return a non-zero status and say so;
  everything else in the library is unaffected. `gfnff_parametrisation_io_available()`
  reports which build you have, `c_gfnff_toml_available()` from C, and
  `gfnff._lib.toml_available()` from Python. Meson leaves toml-f `auto`; CMake
  defaults `WITH_TOMLF` to `ON`, because that build produces the shared library
  the Python bindings load and a wheel without the reader could not be fixed
  from the calling side.

## Choosing a version

The parameter file overlays whichever internal set the *version* selected, so
the two arguments compose. Versions:

| name | effect |
| --- | --- |
| `angewChem2020`, `angewChem2020_1`, `angewChem2020_2` | identical in practice — see below |
| `harmonic2020` | harmonic bonds, for 2D→3D conversion; returns before the EEQ solve, so it has no charges |
| `mcgfnff2023` | molecular crystals; **periodic systems only** |
| `conformer2020` | `angewChem2020_2` with a bond term that cannot dissociate |

`conformer2020` shares the `angewChem2020` parameter set outright — it changes
the *shape* of the bond well past its inflection point, not any parameter — so
a file written for `angewChem2020_2` overlays it unchanged.

The first three are separate enum entries but nothing branches on them, here or
in xtb: the amide and topological-charge fixes their comments describe were
folded in unconditionally. They are kept because published input files name
them, and because a future parametrisation may need the distinction back.
