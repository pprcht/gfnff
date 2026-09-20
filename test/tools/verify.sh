#!/usr/bin/env bash
# This file is part of gfnff.
# SPDX-Identifier: LGPL-3.0-or-later
#
# Refactor gate. Builds, runs the unit tests, then produces the full precision
# reference dump and diffs it against the recorded gold.
#
#   test/tools/verify.sh              build, test, dump, diff against gold
#   test/tools/verify.sh --accept     same, but record the current dump as gold
#   test/tools/verify.sh --tests-only skip the dump, just build and test
#
# The unit tests assert tolerances. The dump asserts bit identity. A structural
# refactor is only clean when both pass; the tests alone would let a change
# through that stays inside every tolerance but still moved the numbers.
set -uo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
BUILD="${GFNFF_BUILD:-$ROOT/build}"
TOOLS="$ROOT/test/tools"
GOLD="$TOOLS/gold/refdump.txt"
SUITES=(gfnff pbc-ini pbc-sp pbc-kernels hessian solvation param-io)

accept=0
tests_only=0
for a in "$@"; do
  case "$a" in
    --accept) accept=1 ;;
    --tests-only) tests_only=1 ;;
    *) echo "unknown option: $a" >&2; exit 2 ;;
  esac
done

if [ ! -d "$BUILD" ]; then
  echo "no build directory at $BUILD (set GFNFF_BUILD to override)" >&2
  exit 2
fi

# ── build lists ───────────────────────────────────────────────────────────────
# cheap, and catches the failure mode where a file was added to one build
# system and not the other; that gap is invisible until something links
echo "== build lists"
"$TOOLS/checkbuild.py" | sed 's/^/   /' || exit 1

# ── build ─────────────────────────────────────────────────────────────────────
echo "== build"
if ! ninja -C "$BUILD" -j"${NJOBS:-8}" >/dev/null 2>&1; then
  # a serial rebuild gives a readable error; the parallel one interleaves
  ninja -C "$BUILD" -j1 2>&1 | tail -40
  echo "BUILD FAILED"
  exit 1
fi
echo "   ok"

# ── unit tests ────────────────────────────────────────────────────────────────
echo "== tests"
total=0
failed=0
for s in "${SUITES[@]}"; do
  out="$("$BUILD/test/gfnff-tester" "$s" 2>&1)"
  rc=$?
  n=$(grep -c '\[PASSED\]' <<<"$out")
  f=$(grep -c '\[FAILED\]' <<<"$out")
  # A suite that aborts (error stop, uncaught error, segfault) stops emitting
  # result lines partway through, so counting them alone reports the tests
  # that ran, not the ones that exist, and a crash reads as a clean pass.
  if [ "$rc" -ne 0 ] && [ "$f" -eq 0 ]; then
    f=1
    out="$out"$'\n''[FAILED] suite aborted with exit code '"$rc"' after '"$n"' tests'
  fi
  total=$((total + n))
  failed=$((failed + f))
  printf '   %-12s %2d passed  %d failed\n' "$s" "$n" "$f"
  [ "$f" -gt 0 ] && grep '\[FAILED\]' <<<"$out" | sed 's/^/      /'
done
echo "   TOTAL $total passed, $failed failed"
[ "$failed" -gt 0 ] && { echo "TESTS FAILED"; exit 1; }
[ "$tests_only" -eq 1 ] && exit 0

# ── reference dump ────────────────────────────────────────────────────────────
DUMP="$BUILD/test/gfnff-refdump"
if [ ! -x "$DUMP" ]; then
  echo "== dump: gfnff-refdump not built, skipping" >&2
  exit 0
fi

echo "== dump"
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
OMP_NUM_THREADS=1 "$DUMP" >"$tmp/t1.txt" || { echo "DUMP FAILED"; exit 1; }
OMP_NUM_THREADS=8 "$DUMP" >"$tmp/t8.txt" || { echo "DUMP FAILED"; exit 1; }

echo "   thread invariance (1 vs 8):"
if "$TOOLS/pdiff.py" "$tmp/t1.txt" "$tmp/t8.txt" --quiet; then
  echo "   ok"
else
  echo "THREAD INVARIANCE BROKEN"
  exit 1
fi

mkdir -p "$(dirname "$GOLD")"
if [ "$accept" -eq 1 ]; then
  cp "$tmp/t1.txt" "$GOLD"
  echo "   recorded new gold: $GOLD"
  exit 0
fi

if [ ! -f "$GOLD" ]; then
  echo "   no gold recorded yet; run with --accept to create one"
  exit 0
fi

echo "   against gold:"
if diff -q "$GOLD" "$tmp/t1.txt" >/dev/null; then
  echo "   bit-identical to gold"
  exit 0
fi
if "$TOOLS/pdiff.py" "$GOLD" "$tmp/t1.txt" --quiet; then
  echo "   within tolerance of gold, but not bit-identical"
  exit 0
fi
echo "DEVIATION FROM GOLD"
exit 1
