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

// Declare the deallocator
extern void c_gfnff_calculator_deallocate(c_gfnff_calculator *calculator);

// Declate the singlepoint calculator
extern void c_gfnff_calculator_singlepoint(c_gfnff_calculator *calculator,
                                              int nat, int *at,
                                              double (*xyz)[3], double *energy,
                                              double (*gradient)[3],
                                              int *iostat);

// Declate the print routine
extern void c_gfnff_calculator_results(c_gfnff_calculator *calculator,
                                       int iunit);


#ifdef __cplusplus
}
#endif

#endif /* GFNFF_INTERFACE_C_H */
