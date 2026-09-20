# Refactor and verification tooling

The unit tests in `test/` assert tolerances. These tools assert bit identity.
Both are needed: a structural change can stay inside every test tolerance and
still have moved the numbers, and only the dump will say so.

## `verify.sh` — the gate

```sh
test/tools/verify.sh              # build, run the suites, dump, diff against gold
test/tools/verify.sh --accept     # same, but record the current dump as the gold
test/tools/verify.sh --tests-only # skip the dump
```

Steps, in order, stopping at the first failure:

1. `ninja` the build directory (`$GFNFF_BUILD`, default `build/`)
2. run all six test suites and report per-suite pass/fail counts
3. run `gfnff-refdump` at 1 and 8 threads and check the two agree — this is
   the OpenMP reduction check that no unit test performs
4. diff the serial dump against `gold/refdump.txt`

Exit status is 0 only when every step passes, so it can be used directly in a
loop or a pre-commit hook.

## `refdump.f90` — the dump

Built as `gfnff-refdump`. Prints energies, gradients, stress tensors and
selected internal arrays at full precision. It asserts nothing, which is why it
is not registered with `test()`.

Three sections:

| section | covers |
|---|---|
| `MOL` | molecular path: hydrogen and halogen bond lists, bond pair matrix, packed non-bonded exponents |
| `PBC` | periodic path: energy, gradient, stress |
| `SOLV` | ALPB internals — Born radii, SASA, derivative matrices |

Systems are walked over three jittered geometries so the list-reuse and
list-capacity paths are exercised, not only the first call. The jitter is a
fixed analytic displacement, never a random number, so two builds see exactly
the same geometries.

Several lines carry an order-sensitive checksum (`chkg`, `chkq`, `brad`, …)
alongside the plain sum. A permuted list or a changed summation order leaves
the sum untouched and moves the checksum, which is the failure mode a
reordering refactor is most likely to introduce.

## `pdiff.py` — the comparison

```sh
test/tools/pdiff.py gold.txt new.txt [--tol 1e-12] [--quiet]
```

Reports the largest deviation per tag, scaled by the largest magnitude under
that tag, so the near-zero net-force entries cannot dominate the ratio. Exits
non-zero if any tag exceeds the tolerance.

## `mkindex.py` — the source index

```sh
test/tools/mkindex.py           # regenerate project_index.md
test/tools/mkindex.py --check   # exit 1 if it is out of date
```

`project_index.md` is generated, not maintained, so it cannot drift from the
tree. Rerun it after anything that moves code between files.

## Recording a new gold

Only when a change is *intended* to move the numbers. Record it in the same
commit as the change, never separately, so that `git log` on the gold file
reads as the list of every deliberate numerical change the project has made.
