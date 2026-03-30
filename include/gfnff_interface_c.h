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

// Declare the deallocator
extern void c_gfnff_calculator_deallocate(c_gfnff_calculator *calculator);

// Declare the singlepoint calculator
// sigma[3][3] receives the stress tensor in Hartree (zero for non-PBC systems)
extern void c_gfnff_calculator_singlepoint(c_gfnff_calculator *calculator,
                                              int nat, int *at,
                                              double (*xyz)[3], double *energy,
                                              double (*gradient)[3],
                                              double sigma[3][3],
                                              int *iostat);

// Declate the print routine
extern void c_gfnff_calculator_results(c_gfnff_calculator *calculator,
                                       int iunit);


#ifdef __cplusplus
}
#endif

#endif /* GFNFF_INTERFACE_C_H */
