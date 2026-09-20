#ifndef GFNFF_INTERFACE_C_H
#define GFNFF_INTERFACE_C_H

#include <stdbool.h> // For booleans
#include <stdio.h>   // For FILE*
#include <stdlib.h>  // For standard library functions


#ifdef __cplusplus
extern "C" {
#endif

// Declare the Fortran structure
typedef struct {
  void *ptr;
} c_gfnff_calculator;

// Declate the initializer
extern c_gfnff_calculator
c_gfnff_calculator_init(int nat, int *at, double (*xyz)[3],
                        int ichrg, int printlevel, const char *solvent);

// Declare the PBC-aware initializer
// lattice[3][3]: three lattice vectors in Bohr (row-major in C)
// npbc: number of periodic dimensions (0-3)
extern c_gfnff_calculator
c_gfnff_calculator_init_pbc(int nat, int *at, double (*xyz)[3],
                             int ichrg, int printlevel,
                             double lattice[3][3], int npbc);

// Parametrisation versions. These mirror TGFFVersionEnum in
// src/core/param_tables.f90; prefer c_gfnff_version_from_name below, which
// asks the library rather than trusting these to stay in step.
#define GFNFF_VERSION_DEFAULT 0
#define GFNFF_VERSION_ANGEWCHEM2020 1
#define GFNFF_VERSION_ANGEWCHEM2020_1 2
#define GFNFF_VERSION_ANGEWCHEM2020_2 3
#define GFNFF_VERSION_MCGFNFF2023 4
#define GFNFF_VERSION_CONFORMER2020 5
#define GFNFF_VERSION_HARMONIC2020 (-1)

// Whether this build can read TOML parameter files (1 = yes, 0 = no).
// Lets a caller tell "built without toml-f" apart from "file missing or
// malformed"; both otherwise surface only as a failed initialization.
extern int c_gfnff_toml_available(void);

// Resolve a version name ("angewChem2020_2", "mcgfnff2023", ...) to its
// integer. Returns GFNFF_VERSION_DEFAULT (0) for an unknown name, which the
// initializer reads as "use the library default".
extern int c_gfnff_version_from_name(const char *name);

// Declare the extended initializer.
// Superset of the two initializers above; every extra input is optional:
//   solvent    : ALPB solvent name; "" or NULL -> gas phase
//   lattice    : double[3][3] lattice vectors (Bohr); NULL -> non-periodic
//   npbc       : number of periodic dimensions (only used with lattice)
//   fraglist   : int[nat] user fragment index per atom; NULL -> automatic
//   refq       : double[nat] atomic reference charges; NULL -> none
//   accuracy   : cutoff/precision factor. Larger is looser and faster; above
//                1.0 the EEQ system is solved in single precision. Pass a
//                value <= 0 for the library default (0.1, or 2.0 above
//                10000 atoms).
//   version    : one of the GFNFF_VERSION_* values above; 0 -> library
//                default. Note that GFNFF_VERSION_MCGFNFF2023 is rejected
//                for non-periodic systems.
//   parametrisation : path to a parameter file; NULL or "" -> the internal
//                set for the selected version. A ".toml" path is applied as
//                an overlay on that set, so the file need only name the keys
//                it changes; any other path is read as the legacy flat
//                format. Requires a build with toml-f for the former.
//   bondmat    : int[nat*nat] molecular graph of integer bond orders; NULL ->
//                GFN-FF perceives the bonds from the geometry. A nonzero
//                element (i,j) declares a bond between atoms i and j, and
//                replaces bond perception entirely, so the connectivity need
//                not be supported by the coordinates. Storage order is
//                irrelevant: the matrix must be symmetric, with a zero
//                diagonal, non-negative entries and at most 41 bonds per
//                atom, all of which are checked. Only the zero/nonzero
//                pattern is used today; the orders are recorded but have no
//                consumer, since GFN-FF derives pi bond orders itself.
//                Molecular systems only (npbc = 0).
extern c_gfnff_calculator
c_gfnff_calculator_init_ex(int nat, int *at, double (*xyz)[3],
                           int ichrg, int printlevel, const char *solvent,
                           double (*lattice)[3], int npbc,
                           int *fraglist, double *refq,
                           double accuracy, int version,
                           const char *parametrisation,
                           const int *bondmat);

// Declare the deallocator
extern void c_gfnff_calculator_deallocate(c_gfnff_calculator *calculator);

// Declare the singlepoint calculator
// sigma[3][3] receives the stress tensor in Hartree (zero for non-PBC systems)
// lattice[3][3]: updated lattice vectors (Bohr); pass NULL to reuse stored lattice
extern void c_gfnff_calculator_singlepoint(c_gfnff_calculator *calculator,
                                              int nat, int *at,
                                              double (*xyz)[3], double *energy,
                                              double (*gradient)[3],
                                              double sigma[3][3],
                                              double (*lattice)[3],
                                              int *iostat);

// Returned in iostat when the Hessian is requested for a periodic system.
// Distinct from the library's own status codes, which are 0 and +/-1.
#define GFNFF_C_PBC_UNSUPPORTED (-2)

// Declare the Hessian
//
// hessian: caller-owned, must already hold 9*nat*nat doubles, and is written
//   only when iostat comes back 0. The matrix is symmetrised before it is
//   returned, so it reads the same row-major here as it does column-major in
//   Fortran -- no transpose is needed. Degree of freedom (c,A) is at
//   3*A + c with zero-based c in [0,3) and zero-based atom index A.
// energy, gradient: pass NULL to skip either.
// step: finite-difference step in Bohr for the terms not yet available in
//   closed form; pass <= 0 to let the library choose.
// iostat: 0 on success, GFNFF_C_PBC_UNSUPPORTED if the calculator was
//   initialised periodic, otherwise the library error code.
//
// Periodic systems are rejected rather than silently finite-differenced: the
// analytic terms have no periodic implementation yet, and the fallback would
// ignore the images without saying so.
extern void c_gfnff_calculator_hessian(c_gfnff_calculator *calculator,
                                       int nat, int *at,
                                       double (*xyz)[3],
                                       double *hessian,
                                       double *energy,
                                       double (*gradient)[3],
                                       double step,
                                       int *iostat);

// Returned in iostat when charges are requested before any singlepoint ran.
#define GFNFF_C_NO_CHARGES (-3)

// Declare the charge accessor
//
// charges: caller-owned, must already hold nat doubles, and is written only
//   when iostat comes back 0. Values are the EEQ partial charges in e.
// iostat: 0 on success, GFNFF_C_NO_CHARGES if no singlepoint has run on this
//   calculator yet, or if nat disagrees with the stored system.
//
// The charges are a by-product of the energy evaluation rather than a
// separate model, so they reflect the geometry of the last singlepoint.
extern void c_gfnff_calculator_charges(c_gfnff_calculator *calculator,
                                       int nat, double *charges,
                                       int *iostat);

// Declate the print routine
extern void c_gfnff_calculator_results(c_gfnff_calculator *calculator,
                                       int iunit);


#ifdef __cplusplus
}
#endif

#endif /* GFNFF_INTERFACE_C_H */
