#!/usr/bin/env python3
# This file is part of gfnff.
# SPDX-Identifier: LGPL-3.0-or-later
"""Compare two refdump outputs.

Reports the largest deviation per tag, scaled by the largest magnitude that
appears under that tag -- so the near-zero net-force entries do not dominate
the ratio and drown out a real difference elsewhere.

Exit status is 0 when every tag agrees to within the tolerance, 1 otherwise,
so this can be used directly as a gate in a shell script.

    pdiff.py gold.txt new.txt [--tol 1e-12] [--quiet]
"""
import sys


def num(t):
    try:
        return float(t)
    except ValueError:
        return None


def load(path):
    """tag -> list of every number appearing on lines with that tag."""
    d = {}
    for line in open(path):
        f = line.split()
        if len(f) < 3:
            continue
        d.setdefault(f[0], []).extend(v for v in map(num, f[1:]) if v is not None)
    return d


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    tol = 1e-12
    if "--tol" in sys.argv:
        tol = float(sys.argv[sys.argv.index("--tol") + 1])
        args = [a for a in args if a != str(tol)]
    quiet = "--quiet" in sys.argv

    if len(args) < 2:
        print(__doc__)
        return 2

    a, b = load(args[0]), load(args[1])
    bad = 0
    worst_overall = 0.0

    for k in sorted(a):
        if k not in b:
            print(f"  {k}: MISSING in {args[1]}")
            bad = 1
            continue
        if len(a[k]) != len(b[k]):
            print(f"  {k}: field count {len(a[k])} vs {len(b[k])}")
            bad = 1
            continue
        scale = max(max(map(abs, a[k]), default=0.0), 1e-12)
        worst = max((abs(x - y) for x, y in zip(a[k], b[k])), default=0.0)
        rel = worst / scale
        worst_overall = max(worst_overall, rel)
        if rel > tol:
            print(f"  {k:14s} max abs dev {worst:.3e}   (scaled {rel:.3e})  <-- DIFFERS")
            bad = 1
        elif not quiet:
            print(f"  {k:14s} max abs dev {worst:.3e}   (scaled {rel:.3e})")

    for k in sorted(b):
        if k not in a:
            print(f"  {k}: EXTRA in {args[1]}")
            bad = 1

    print(f"  worst scaled deviation over all tags: {worst_overall:.3e}  (tol {tol:.1e})")
    return bad


if __name__ == "__main__":
    sys.exit(main())
