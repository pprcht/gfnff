# This file is part of gfnff.
# SPDX-Identifier: LGPL-3.0-or-later
"""
Command-line interface for GFN-FF singlepoint and geometry optimisation.

Usage
-----
  gfnff <input> [options]

The input file is read by ASE (any supported format: xyz, extxyz, POSCAR, cif, …).
"""

import argparse
import sys

import numpy as np


_HARTREE_TO_EV = 27.211386245988   # matches ase.units.Hartree
_BOHR_TO_ANG = 0.529177210903       # matches ase.units.Bohr


def _parse_args(argv=None):
    parser = argparse.ArgumentParser(
        prog="gfnff",
        description=(
            "GFN-FF singlepoint energy / geometry optimisation.\n"
            "The input structure is read by ASE (xyz, extxyz, POSCAR, cif, …)."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "input",
        help="Input structure file.",
    )
    parser.add_argument(
        "--chrg",
        type=int,
        default=0,
        metavar="INT",
        help="Total molecular charge (default: 0).",
    )

    # Optimisation
    opt_group = parser.add_argument_group("optimisation")
    opt_group.add_argument(
        "--opt",
        action="store_true",
        help=(
            "Perform geometry optimisation (L-BFGS). "
            "Writes the optimisation path to gfnff.log.extxyz."
        ),
    )
    opt_group.add_argument(
        "--fmax",
        type=float,
        default=0.01,
        metavar="FLOAT",
        help="Maximum force component for convergence in eV/Å (default: 0.01).",
    )
    opt_group.add_argument(
        "--maxsteps",
        type=int,
        default=500,
        metavar="INT",
        help="Maximum number of optimisation steps (default: 500).",
    )
    opt_group.add_argument(
        "--outfile",
        default="gfnff.log.extxyz",
        metavar="FILE",
        help="Optimisation trajectory output file (default: gfnff.log.extxyz).",
    )

    # Solvation: --alpb and --solv are aliases
    solv_group = parser.add_mutually_exclusive_group()
    solv_group.add_argument(
        "--alpb",
        metavar="SOLVENT",
        help=(
            "Activate implicit solvation and select the solvent "
            "(e.g. h2o, acetone, chcl3, dmso, thf, …)."
        ),
    )
    solv_group.add_argument(
        "--solv",
        metavar="SOLVENT",
        help="Alias for --alpb.",
    )

    return parser.parse_args(argv)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _print_header(atoms, solvent, charge):
    nat = len(atoms)
    formula = atoms.get_chemical_formula()
    pbc_str = (
        "3D periodic" if all(atoms.pbc)
        else f"PBC {atoms.pbc.tolist()}" if any(atoms.pbc)
        else "molecular"
    )
    print("=" * 60)
    print("  GFN-FF  —  standalone Python interface")
    print("=" * 60)
    print(f"  Formula     : {formula}")
    print(f"  Atoms       : {nat}")
    print(f"  Charge      : {charge:+d}")
    print(f"  System type : {pbc_str}")
    if solvent:
        print(f"  Solvent     : {solvent} (ALPB/GBSA)")
    print("=" * 60)
    print()


def _print_singlepoint_results(atoms):
    from ase.units import Bohr, Hartree

    energy_ev = atoms.get_potential_energy()
    energy_eh = energy_ev / Hartree
    forces = atoms.get_forces()             # eV/Å
    # convert to Eh/Bohr for the gnorm printout
    grad_eh_bohr = -forces * (Bohr / Hartree)
    gnorm = np.linalg.norm(grad_eh_bohr, axis=1).mean()

    print(f"  Total energy : {energy_eh:18.10f} Eh")
    print(f"               : {energy_ev:18.10f} eV")
    print(f"  Gradient norm: {gnorm:18.10f} Eh/a0  (mean |∇E|)")
    print()


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def main(argv=None):
    args = _parse_args(argv)
    solvent = args.alpb or args.solv or ""

    # --- ASE import check -----------------------------------------------
    try:
        import ase.io
    except ImportError:
        print(
            "Error: ASE is required for the gfnff CLI.\n"
            "Install it with:  pip install 'gfnff[ase]'",
            file=sys.stderr,
        )
        sys.exit(1)

    from gfnff import GFNFF

    # --- Read structure --------------------------------------------------
    try:
        atoms = ase.io.read(args.input)
    except Exception as exc:
        print(f"Error reading '{args.input}': {exc}", file=sys.stderr)
        sys.exit(1)

    atoms.info["charge"] = args.chrg

    calc = GFNFF(charge=args.chrg, solvent=solvent, printlevel=1)
    atoms.calc = calc

    _print_header(atoms, solvent, args.chrg)

    # --- Singlepoint or optimisation ------------------------------------
    if args.opt:
        _run_opt(atoms, args)
    else:
        _run_singlepoint(atoms)


def _run_singlepoint(atoms):
    print("  Singlepoint calculation")
    print("-" * 60)
    try:
        _print_singlepoint_results(atoms)
        atoms.calc._gfnff.print_results()
    except Exception as exc:
        print(f"Calculation failed: {exc}", file=sys.stderr)
        sys.exit(1)


def _run_opt(atoms, args):
    import ase.io
    from ase.optimize import LBFGS

    outfile = args.outfile
    print(f"  Geometry optimisation  (fmax = {args.fmax} eV/Å)")
    print(f"  Trajectory → {outfile}")
    print("-" * 60)

    # Write each accepted step to extxyz; the first call truncates the file.
    step = [0]

    def _write_frame():
        # Copy energy and forces into atoms.info / atoms.arrays so the
        # extxyz writer picks them up automatically.
        res = atoms.calc.results
        frame = atoms.copy()
        frame.info["energy"] = res["energy"]
        frame.arrays["forces"] = res["forces"].copy()
        # Append after the first frame so the file is not truncated mid-run.
        ase.io.write(outfile, frame, format="extxyz", append=(step[0] > 0))
        step[0] += 1

    opt = LBFGS(atoms, logfile="-")
    opt.attach(_write_frame, interval=1)

    try:
        converged = opt.run(fmax=args.fmax, steps=args.maxsteps)
    except Exception as exc:
        print(f"\nOptimisation failed: {exc}", file=sys.stderr)
        sys.exit(1)

    print()
    if converged:
        print("  Converged.")
    else:
        print(f"  WARNING: optimisation did not converge in {args.maxsteps} steps.")

    print()
    _print_singlepoint_results(atoms)
    print(f"  Trajectory written to '{outfile}'  ({step[0]} frames)")
    print()
