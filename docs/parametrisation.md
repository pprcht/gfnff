# Force-field versions and molecular graphs

Everything on this page applies to both Python calculators, `GFNFFCalculator` and the
ASE calculator `GFNFF` (see [python.md](python.md)).

## Choosing a parametrisation

Both calculators take `version` and `parametrisation`.

```python
from gfnff import GFNFFCalculator, Version

# a different force-field version
calc = GFNFFCalculator(numbers, positions, version="harmonic2020")
calc = GFNFFCalculator(numbers, positions, version=Version.harmonic2020)  # same

# a custom parameter file, overlaid on the internal set for that version
calc = GFNFFCalculator(numbers, positions, parametrisation="my-set.toml")
```

| `version` | Effect |
|-----------|--------|
| `None` *(default)* | `angewChem2020_2` |
| `angewChem2020`, `angewChem2020_1`, `angewChem2020_2` | identical in practice; nothing branches on the distinction, here or in xtb |
| `harmonic2020` | harmonic bond potential for 2D→3D conversion; runs no EEQ solve, so `charges()` raises |
| `mcgfnff2023` | molecular crystals; rejected for non-periodic systems |
| `conformer2020` | `angewChem2020_2` with bonds that cannot dissociate, see [below](#conformer2020-bonds-that-cannot-break) |

A `.toml` parameter file is an **overlay**: it need only name the keys it
changes, and the rest keep the internal values for the selected version. See
[`param/README.md`](../param/README.md) for the format and for how to dump the
current set as a starting point. TOML support needs a build with toml-f
(`gfnff._lib.toml_available()`).

### `conformer2020`: bonds that cannot break

The published bond term is a Gaussian well in the deviation from a
CN-dependent reference length. It has an inflection point at about 0.5 Å of
stretch and decays to zero beyond it, so a bond can be pulled apart for a
finite cost and the molecule can dissociate during an optimisation or an MD
run while the bond list still contains the bond.

`conformer2020` keeps that well exactly where it is convex and continues it
past the inflection point along its own tangent, which rises without bound at
a constant restoring force. The join is C², so gradients and Hessians stay
continuous.

Inside the convex region the results are bit-identical to `angewChem2020_2`,
since the same code path runs. At the GFN-FF minimum of caffeine
the most stretched bond sits 0.17 Å from its reference, against a switch at
0.49 Å, so ordinary conformers, thermal MD and normal strain never reach the
continuation:

```python
atoms.calc = GFNFF(version="conformer2020")   # same energies, no dissociation
```

Pulling one C–H bond out with the topology held fixed, energies relative to
the minimum in eV:

| stretch | `angewChem2020_2` | `conformer2020` |
|--------:|------------------:|----------------:|
| 0.2 Å | 0.54 | 0.54 |
| 0.4 Å | 1.93 | 1.93 |
| 1.0 Å | 4.18 | 5.33 |
| 4.0 Å | 4.52 | 20.77 |
| 8.0 Å | 4.52 | 41.35 |

The published term saturates at its well depth; the variant climbs at
5.15 eV/Å, which is `sqrt(2α)·D·exp(-1/2)` for that bond.

All other terms are unchanged. The angle and torsion terms are damped to zero
as their bonds stretch, but the damping is never reached once the bonds cannot
stretch that far, so these terms need no modification.

## Supplying the molecular graph

By default GFN-FF works out the connectivity from the geometry. `bond_matrix`
replaces that step with a graph you supply, an `(nat, nat)` integer matrix
where a nonzero element declares a bond:

```python
calc = GFNFFCalculator(numbers, positions, bond_matrix=bm)
atoms.calc = GFNFF(bond_matrix=bm)          # or atoms.info["bond_matrix"] = bm
```

Handing back the graph GFN-FF would have perceived reproduces the same force
field exactly, i.e. `bond_matrix` overrides the perception and is not a
separate mode. It works with every version.

The matrix must be symmetric, zero on the diagonal, non-negative, and carry at
most 41 bonds per atom. All of this is checked, and a violation is an error.
Molecular systems only, since a graph has no cell index. Only the zero/nonzero
pattern is currently used. The magnitudes are validated and stored, but not
read: GFN-FF derives its own π bond orders from a Hückel
treatment, and `harmonic2020` uses the bond list alone.

### Building a structure from a graph

The option exists because with `harmonic2020` the entire force field is
determined by the graph and the elements. Its bond term targets
`0.7·(rcov_i + rcov_j)`, its topological charges come from a shortest-path walk
over the graph, and its pair exponents follow from those. The geometry enters
only through the repulsion. Random atom positions plus a graph are therefore
sufficient to build a structure:

```python
from ase.optimize import BFGS
from gfnff import GFNFF

soup.calc = GFNFF(version="harmonic2020", bond_matrix=bm)   # random positions
BFGS(soup).run(fmax=1e-3)                                   # -> 3D structure
soup.calc = GFNFF()                                         # then regular GFN-FF
BFGS(soup).run(fmax=1e-3)
```

Starting from uniformly random coordinates, caffeine is recovered with the correct
bond lengths (1.261 ± 0.205 Å against 1.261 ± 0.152 Å for the reference), and
relaxing that with ordinary GFN-FF reaches the same minimum to within 0.002 eV.

When a graph is supplied *and* the version is `harmonic2020`, the perception
phases are skipped entirely: rings, hybridisation, π systems and the bonded
parameters all read the geometry, which in this mode carries no information.

