#include "gfnff_interface_c.h"

// The test implementation program
int main() {
  // Test molecule: caffeine
  int nat = 24; // Number of atoms
  int at[24] = {6, 7, 6, 7, 6, 6, 6, 8, 7, 6, 8, 7,
                6, 6, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}; // Atom types
  double xyz[24][3] = {
      {2.02799738646442, 0.09231312124713, -0.14310895950963},
      {4.75011007621000, 0.02373496014051, -0.14324124033844},
      {6.33434307654413, 2.07098865582721, -0.14235306905930},
      {8.72860718071825, 1.38002919517619, -0.14265542523943},
      {8.65318821103610, -1.19324866489847, -0.14231527453678},
      {6.23857175648671, -2.08353643730276, -0.14218299370797},
      {5.63266886875962, -4.69950321056008, -0.13940509630299},
      {3.44931709749015, -5.48092386085491, -0.14318454855466},
      {7.77508917214346, -6.24427872938674, -0.13107140408805},
      {10.30229550927022, -5.39739796609292, -0.13672168520430},
      {12.07410272485492, -6.91573621641911, -0.13666499342053},
      {10.70038521493902, -2.79078533715849, -0.14148379504141},
      {13.24597858727017, -1.76969072232377, -0.14218299370797},
      {7.40891694074004, -8.95905928176407, -0.11636933482904},
      {1.38702118184179, 2.05575746325296, -0.14178615122154},
      {1.34622199478497, -0.86356704498496, 1.55590600570783},
      {1.34624089204623, -0.86133716815647, -1.84340893849267},
      {5.65596919189118, 4.00172183859480, -0.14131371969009},
      {14.67430918222276, -3.26230980007732, -0.14344911021228},
      {13.50897177220290, -0.60815166181684, 1.54898960808727},
      {13.50780014200488, -0.60614855212345, -1.83214617078268},
      {5.41408424778406, -9.49239668625902, -0.11022772492007},
      {8.31919801555568, -9.74947502841788, 1.56539243085954},
      {8.31511620712388, -9.76854236502758, -1.79108242206824}};

  // Test variables
  int ichrg = 0;
  int printlevel = 1;
  char *solvent="";

  // Call the Fortran function
  c_gfnff_calculator calc =
      c_gfnff_calculator_init(nat, at, xyz, ichrg, printlevel, solvent);

  if (calc.ptr == NULL) {
    printf("Error initializing gfnff calculator.\n");
    return 1;
  }

  printf("gfnff calculator initialized successfully.\n");

  // Use the calculator...
  double energy;
  double gradient[nat][3]; // Adjust the size to match nat
  int iostat;

  // Call the singlepoint function
  c_gfnff_calculator_singlepoint(&calc, nat, at, xyz, &energy, gradient,
                                 &iostat);

  // Check the result and print it
  if (iostat == 0) {
    printf("Singlepoint calculation successful.\n");
    printf("Energy: %f\n", energy);

    // Print the gradient for the first atoms(optional)
    for (int i = 0; i < 3; i++) {
      int j = 0;
      printf("Gradient[%d][%d] = %e\n", j, i, gradient[j][i]);
    }
  } else {
    printf("Singlepoint calculation failed with iostat = %d\n", iostat);
  }

  // Also test the intrinsic print routine
  int iunit = 6; // Use 6 to get Fortran's STDOUT

  // Call the print function
  c_gfnff_calculator_results(&calc, iunit);

  // Deallocate the Fortran object
  c_gfnff_calculator_deallocate(&calc);

  if (calc.ptr == NULL) {
    printf("gfnff calculator deallocated successfully.\n");
  } else {
    printf("Error deallocating gfnff calculator.\n");
  }

  return 0;
}
