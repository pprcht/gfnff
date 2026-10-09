#!/usr/bin/env python3
# This file is part of gfnff.
# SPDX-Identifier: LGPL-3.0-or-later
"""Assert that the meson and CMake source lists agree with each other and with
the tree.

Both build systems list every source by hand, and they have drifted apart
before -- a file was added to one and not the other, and the gap went unnoticed
because each build works fine on its own until someone links against the
missing symbol. This check makes that drift a hard error.

Three things are verified:

  1. every ``.f90``/``.F90`` under ``src/`` appears in the meson lists
  2. the same holds for the CMake lists
  3. nothing is listed that does not exist on disk

    test/tools/checkbuild.py
"""
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
SRC = os.path.join(ROOT, "src")

RE_MESON = re.compile(r"'([^']+\.[fF]90)'")
RE_CMAKE = re.compile(r'"\$\{dir\}/([^"]+\.[fF]90)"')


def on_disk():
    found = set()
    for dirpath, _, filenames in os.walk(SRC):
        for fn in filenames:
            if fn.endswith((".f90", ".F90")):
                found.add(os.path.relpath(os.path.join(dirpath, fn), SRC))
    return found


def listed(buildfile, pattern):
    """Collect sources named by every build file of the given name under src/."""
    found = set()
    for dirpath, _, filenames in os.walk(SRC):
        if buildfile not in filenames:
            continue
        d = os.path.relpath(dirpath, SRC)
        with open(os.path.join(dirpath, buildfile)) as fh:
            for name in pattern.findall(fh.read()):
                found.add(os.path.normpath(os.path.join(d, name)))
    return found


def report(label, missing, extra):
    bad = False
    for f in sorted(missing):
        print(f"  {label}: {f} exists on disk but is not listed")
        bad = True
    for f in sorted(extra):
        print(f"  {label}: {f} is listed but does not exist")
        bad = True
    return bad


def main():
    disk = on_disk()
    meson = listed("meson.build", RE_MESON)
    cmake = listed("CMakeLists.txt", RE_CMAKE)

    bad = False
    bad |= report("meson", disk - meson, meson - disk)
    bad |= report("cmake", disk - cmake, cmake - disk)

    only_meson = meson - cmake
    only_cmake = cmake - meson
    for f in sorted(only_meson):
        print(f"  drift: {f} is in meson but not in CMake")
        bad = True
    for f in sorted(only_cmake):
        print(f"  drift: {f} is in CMake but not in meson")
        bad = True

    if bad:
        print("BUILD LISTS INCONSISTENT")
        return 1
    print(f"build lists agree: {len(disk)} sources in meson, CMake and on disk")
    return 0


if __name__ == "__main__":
    sys.exit(main())
