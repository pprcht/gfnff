# Performance

![GFN-FF benchmark](../assets/benchmark.png)

Caffeine clusters from 24 to 1536 atoms, gfortran 14.3 `-O3`, pinned to 8
physical cores, OpenBLAS 0.3.31 restricted to one thread so that OpenMP provides
the parallelism, best of 3 runs. Both trees give the same energies to 2e-16.

- **Energy and gradient** are 1.2–1.6x faster than the pre-refactor tree
  from 192 atoms upwards. Below that a call takes about a millisecond or less
  and the ratio is dominated by timer noise.
- **The analytic Hessian** is 28x faster than a full 6N-gradient finite
  difference at 24 atoms and 59x faster at 768 atoms (1.8 s instead of 107 s).

The setup is unchanged by the refactor and remains the practical size limit: at
1536 atoms topology perception takes ~2.9 s against ~102 ms for a gradient. The
energy and gradient speedup holds on a single thread as well: at 1536 atoms the
current tree takes 174 ms against 248 ms, and it scales to 1.7x on 8 cores where
the old tree reaches 1.5x.

The analytic Hessian is dominated by three large `DGEMM` calls and one
factorise-and-solve step. The electrostatic, dispersion and bond terms together
account for ~93% of it, and each contracts a `(3N,N)` derivative block into the
`(3N,3N)` Hessian. Its OpenMP parallelisation therefore gains only ~1.2x on 8
cores, whereas a threaded BLAS gives ~3.5x. Hessian performance hence depends
mainly on the linear algebra backend. Reference netlib LAPACK is serial and
costs roughly 1.6x on singlepoints and considerably more on Hessians.
`find_package(LAPACK)` takes the first library it finds; pass
`-DBLA_VENDOR=...` (e.g. `OpenBLAS`, `Intel10_64lp`) to select one explicitly.

The benchmark harness is in `bench/` (untracked, it generates its own
structures); `bench/build.sh` builds both trees, `bench/run2.sh` reproduces the
raw numbers and `bench/plot.py` redraws the figure.
