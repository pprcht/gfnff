#include <math.h>
#include <stdlib.h>
#include <stdio.h>
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
  double sigma[3][3];      // stress tensor (zero for non-PBC)
  int iostat;

  // Call the singlepoint function (NULL lattice: non-PBC, reuse stored)
  c_gfnff_calculator_singlepoint(&calc, nat, at, xyz, &energy, gradient,
                                 sigma, NULL, &iostat);

  // Check the result and print it
  if (iostat == 0) {
    printf("Singlepoint calculation successful.\n");
    printf("Energy: %f\n", energy);

    // Print the gradient for the first atoms(optional)
    for (int i = 0; i < 3; i++) {
      int j = 0;
      printf("Gradient[%d][%d] = %e\n", j, i, gradient[j][i]);
    }

    // Print the stress tensor
    printf("Sigma (molecular, should be zeroed):\n");
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        printf("  sigma[%d][%d] = %e\n", i, j, sigma[i][j]);
  } else {
    printf("Singlepoint calculation failed with iostat = %d\n", iostat);
  }

  // ── charges ─────────────────────────────────────────────────────────────
  // Caller-owned buffer, written only on iostat 0. The EEQ charges must sum
  // to the total charge given at init, which is the one property of theirs
  // that holds for any geometry and so can be asserted without a reference.
  {
    double charges[nat];
    int q_iostat;
    c_gfnff_calculator_charges(&calc, nat, charges, &q_iostat);
    if (q_iostat != 0) {
      printf("Charge lookup failed with iostat = %d\n", q_iostat);
      return 1;
    }
    double qsum = 0.0;
    for (int i = 0; i < nat; i++)
      qsum += charges[i];
    printf("Charges: q[0] = %f, sum = %e (expected %d)\n", charges[0], qsum,
           ichrg);
    if (fabs(qsum - (double)ichrg) > 1.0e-8) {
      printf("FAIL: charges do not sum to the total charge\n");
      return 1;
    }
  }

  // Also test the intrinsic print routine
  int iunit = 6; // Use 6 to get Fortran's STDOUT

  // Call the print function
  c_gfnff_calculator_results(&calc, iunit);

  // Deallocate the Fortran object
  // ── Hessian ─────────────────────────────────────────────────────────────
  // 3N x 3N, caller-owned. Checked two ways: the matrix must be symmetric
  // (the library symmetrises, so any asymmetry here means the row-major
  // handoff is wrong), and one column must match a central difference of
  // the gradient the C API itself returns.
  {
    const int n3 = 3 * nat;
    double *hess = malloc((size_t)n3 * n3 * sizeof(double));
    if (hess == NULL) {
      printf("Hessian allocation failed.\n");
      return 1;
    }
    double h_energy;
    double h_grad[24][3];
    int h_iostat;

    c_gfnff_calculator_hessian(&calc, nat, at, xyz, hess, &h_energy, h_grad,
                               0.0, &h_iostat);
    if (h_iostat != 0) {
      printf("Hessian failed with iostat = %d\n", h_iostat);
      free(hess);
      return 1;
    }
    printf("Hessian computed successfully.\n");
    printf("Hessian energy: %f\n", h_energy);
    printf("Hessian[0][0] = %e\n", hess[0]);
    printf("Hessian[0][1] = %e\n", hess[1]);

    // symmetry
    double asym = 0.0;
    for (int i = 0; i < n3; i++)
      for (int j = 0; j < n3; j++) {
        double d = fabs(hess[(size_t)i * n3 + j] - hess[(size_t)j * n3 + i]);
        if (d > asym) asym = d;
      }
    printf("Hessian max asymmetry: %e\n", asym);
    if (asym > 1e-10) {
      printf("Hessian is not symmetric -- row/column handoff is wrong.\n");
      free(hess);
      return 1;
    }

    // column 0 against a central difference of the API gradient
    const double delta = 5.0e-4;
    double gplus[24][3], gminus[24][3], e_tmp, sig_tmp[3][3];
    int st;
    xyz[0][0] += delta;
    c_gfnff_calculator_singlepoint(&calc, nat, at, xyz, &e_tmp, gplus, sig_tmp,
                                   NULL, &st);
    xyz[0][0] -= 2.0 * delta;
    c_gfnff_calculator_singlepoint(&calc, nat, at, xyz, &e_tmp, gminus, sig_tmp,
                                   NULL, &st);
    xyz[0][0] += delta;

    double worst = 0.0;
    for (int a = 0; a < nat; a++)
      for (int c = 0; c < 3; c++) {
        double fd = (gplus[a][c] - gminus[a][c]) / (2.0 * delta);
        double d = fabs(fd - hess[(size_t)(3 * a + c) * n3 + 0]);
        if (d > worst) worst = d;
      }
    printf("Hessian column 0 vs finite difference, max dev: %e\n", worst);
    if (worst > 1e-4) {
      printf("Hessian disagrees with the finite-differenced gradient.\n");
      free(hess);
      return 1;
    }

    // NULL energy and gradient must be accepted
    c_gfnff_calculator_hessian(&calc, nat, at, xyz, hess, NULL, NULL, 0.0,
                               &h_iostat);
    if (h_iostat != 0) {
      printf("Hessian with NULL outputs failed with iostat = %d\n", h_iostat);
      free(hess);
      return 1;
    }
    printf("Hessian accepts NULL energy and gradient.\n");
    free(hess);
  }

  c_gfnff_calculator_deallocate(&calc);

  if (calc.ptr == NULL) {
    printf("gfnff calculator deallocated successfully.\n");
  } else {
    printf("Error deallocating gfnff calculator.\n");
  }

  // ── version selection ───────────────────────────────────────────────────────
  // The names must resolve to the documented integers, and harmonic2020 must
  // actually reach the force field: it replaces the bond potential outright,
  // so its energy cannot match the default parametrisation.
  {
    printf("\n--- version selection ---\n");
    if (c_gfnff_version_from_name("angewChem2020_2") !=
        GFNFF_VERSION_ANGEWCHEM2020_2) {
      printf("FAIL: angewChem2020_2 does not resolve to its header value\n");
      return 1;
    }
    if (c_gfnff_version_from_name("mcgfnff2023") != GFNFF_VERSION_MCGFNFF2023) {
      printf("FAIL: mcgfnff2023 does not resolve to its header value\n");
      return 1;
    }
    if (c_gfnff_version_from_name("no-such-version") != GFNFF_VERSION_DEFAULT) {
      printf("FAIL: unknown version name does not fall back to the default\n");
      return 1;
    }

    c_gfnff_calculator hcalc = c_gfnff_calculator_init_ex(
        nat, at, xyz, ichrg, 0, "", NULL, 0, NULL, NULL, -1.0,
        c_gfnff_version_from_name("harmonic2020"), NULL, NULL);
    if (hcalc.ptr == NULL) {
      printf("FAIL: could not initialise the harmonic2020 version\n");
      return 1;
    }
    double h_e;
    double h_g[nat][3];
    double h_s[3][3];
    int h_io;
    c_gfnff_calculator_singlepoint(&hcalc, nat, at, xyz, &h_e, h_g, h_s, NULL,
                                   &h_io);
    if (h_io != 0) {
      printf("FAIL: harmonic2020 singlepoint failed with iostat = %d\n", h_io);
      return 1;
    }
    printf("harmonic2020 energy: %f (default was %f)\n", h_e, energy);
    if (fabs(h_e - energy) < 1.0e-6) {
      printf("FAIL: harmonic2020 gave the default energy, version ignored\n");
      return 1;
    }

    // harmonic2020 returns before the EEQ solve, so it must report that it has
    // no charges rather than hand back the zeros left in the buffer
    double h_q[nat];
    int hq_io;
    c_gfnff_calculator_charges(&hcalc, nat, h_q, &hq_io);
    if (hq_io != GFNFF_C_NO_CHARGES) {
      printf("FAIL: harmonic2020 reported charges (iostat = %d)\n", hq_io);
      return 1;
    }
    printf("harmonic2020 correctly reports no charges.\n");
    c_gfnff_calculator_deallocate(&hcalc);

    // mcgfnff2023 is for molecular crystals and must be refused here
    c_gfnff_calculator mcalc = c_gfnff_calculator_init_ex(
        nat, at, xyz, ichrg, 0, "", NULL, 0, NULL, NULL, -1.0,
        GFNFF_VERSION_MCGFNFF2023, NULL, NULL);
    if (mcalc.ptr != NULL) {
      printf("FAIL: mcgfnff2023 was accepted for a non-periodic system\n");
      return 1;
    }
    printf("mcgfnff2023 correctly refused for a non-periodic system.\n");

    // conformer2020 differs from the default only in the dissociative tail of
    // the bond term, so on an undistorted structure it must give back exactly
    // the default energy -- the one version whose selection is invisible here.
    if (c_gfnff_version_from_name("conformer2020") !=
        GFNFF_VERSION_CONFORMER2020) {
      printf("FAIL: conformer2020 does not resolve to its header value\n");
      return 1;
    }
    c_gfnff_calculator ccalc = c_gfnff_calculator_init_ex(
        nat, at, xyz, ichrg, 0, "", NULL, 0, NULL, NULL, -1.0,
        GFNFF_VERSION_CONFORMER2020, NULL, NULL);
    if (ccalc.ptr == NULL) {
      printf("FAIL: could not initialise the conformer2020 version\n");
      return 1;
    }
    double c_e;
    double c_g[nat][3];
    double c_s[3][3];
    int c_io;
    c_gfnff_calculator_singlepoint(&ccalc, nat, at, xyz, &c_e, c_g, c_s, NULL,
                                   &c_io);
    if (c_io != 0) {
      printf("FAIL: conformer2020 singlepoint failed with iostat = %d\n", c_io);
      return 1;
    }
    // The comparison has to be against a calculator built the same way. The
    // one at the top of this file came from c_gfnff_calculator_init, whose
    // defaults differ from init_ex's, and that alone moves the last bits.
    c_gfnff_calculator dcalc = c_gfnff_calculator_init_ex(
        nat, at, xyz, ichrg, 0, "", NULL, 0, NULL, NULL, -1.0,
        GFNFF_VERSION_DEFAULT, NULL, NULL);
    if (dcalc.ptr == NULL) {
      printf("FAIL: could not initialise the default version via init_ex\n");
      return 1;
    }
    double d_e;
    double d_g[nat][3];
    double d_s[3][3];
    int d_io;
    c_gfnff_calculator_singlepoint(&dcalc, nat, at, xyz, &d_e, d_g, d_s, NULL,
                                   &d_io);
    if (d_io != 0) {
      printf("FAIL: default singlepoint failed with iostat = %d\n", d_io);
      return 1;
    }
    // Single threaded the two are bit-identical, but the OpenMP reductions
    // sum in whatever order the threads finish, which moves either energy by
    // a couple of ULP run to run. 1e-12 is far below the ~0.4 Eh the
    // continuation branch is worth once it engages, and far above that noise.
    if (fabs(c_e - d_e) > 1.0e-12) {
      printf("FAIL: conformer2020 moved the energy at equilibrium: %.17g vs "
             "%.17g\n",
             c_e, d_e);
      return 1;
    }
    printf("conformer2020 reproduces the default energy exactly.\n");
    c_gfnff_calculator_deallocate(&dcalc);
    c_gfnff_calculator_deallocate(&ccalc);

    // ── host-supplied molecular graph ────────────────────────────────────
    // A caffeine graph built from the reference geometry, handed back as a
    // bond matrix. Storage order is irrelevant here because the matrix is
    // symmetric, which is also what the library requires.
    {
      int *bm = calloc((size_t)nat * nat, sizeof(int));
      if (bm == NULL) {
        printf("bond matrix allocation failed\n");
        return 1;
      }
      int nbond = 0;
      for (int i = 0; i < nat; i++)
        for (int j = 0; j < i; j++) {
          double d = 0.0;
          for (int c = 0; c < 3; c++) {
            double dx = xyz[i][c] - xyz[j][c];
            d += dx * dx;
          }
          // 3.2 Bohr covers every bond in caffeine and no non-bonded pair
          if (d < 3.2 * 3.2) {
            bm[(size_t)i * nat + j] = 1;
            bm[(size_t)j * nat + i] = 1;
            nbond++;
          }
        }
      printf("supplied graph: %d bonds\n", nbond);

      c_gfnff_calculator gcalc = c_gfnff_calculator_init_ex(
          nat, at, xyz, ichrg, 0, "", NULL, 0, NULL, NULL, -1.0,
          GFNFF_VERSION_DEFAULT, NULL, bm);
      if (gcalc.ptr == NULL) {
        printf("FAIL: could not initialise from a supplied molecular graph\n");
        free(bm);
        return 1;
      }
      double g_e;
      double g_g[nat][3];
      double g_s[3][3];
      int g_io;
      c_gfnff_calculator_singlepoint(&gcalc, nat, at, xyz, &g_e, g_g, g_s, NULL,
                                     &g_io);
      if (g_io != 0) {
        printf("FAIL: graph singlepoint failed with iostat = %d\n", g_io);
        free(bm);
        return 1;
      }
      // The graph matches what perception would find on this geometry, so the
      // energy must too; a mismatch means the matrix was misread.
      if (fabs(g_e - d_e) > 1.0e-12) {
        printf("FAIL: supplied graph moved the energy: %.17g vs %.17g\n", g_e,
               d_e);
        free(bm);
        return 1;
      }
      printf("supplied graph reproduces the perceived force field.\n");
      c_gfnff_calculator_deallocate(&gcalc);

      // an asymmetric graph must be refused, not quietly repaired
      bm[0 * nat + 1] = 1;
      bm[1 * nat + 0] = 0;
      c_gfnff_calculator bad = c_gfnff_calculator_init_ex(
          nat, at, xyz, ichrg, 0, "", NULL, 0, NULL, NULL, -1.0,
          GFNFF_VERSION_DEFAULT, NULL, bm);
      if (bad.ptr != NULL) {
        printf("FAIL: an asymmetric molecular graph was accepted\n");
        free(bm);
        return 1;
      }
      printf("asymmetric graph correctly refused.\n");
      free(bm);
    }
  }

  // ── PBC test: SiO2 alpha-quartz unit cell ────────────────────────────────────
  printf("\n--- PBC test: SiO2 alpha-quartz unit cell ---\n");

  int nat_pbc = 9;
  int at_pbc[9] = {8, 8, 8, 8, 8, 8, 14, 14, 14};
  double xyz_pbc[9][3] = {
      { 2.82781861325240,  2.96439280874170,  3.12827803849279},
      { 7.19124230791576,  0.98723342603994,  4.89004701836746},
      { 4.95491880597601,  4.82830910314898,  8.74847811174740},
      { 0.19290883043307,  2.30645007856310,  8.72969832061507},
      {-2.01592208020090,  6.16478744235115,  4.87273962147340},
      { 0.66183062221384,  7.07392578563696,  0.27767968372345},
      { 4.55701736204879,  0.06291337111965,  3.31745840478609},
      {-2.10064209975148,  3.63969476409878,  6.81014625000326},
      { 2.31009832827224,  4.12572862149043,  0.08842485276656}};
  // Hexagonal lattice: a=b=9.284 Bohr, c=10.214 Bohr, gamma=120 deg
  // Each C row maps to a Fortran column (lattice vector)
  double a_sio2 = 9.28422449595511046;
  double c_sio2 = 10.21434769907115;
  double lattice_sio2[3][3] = {
      {a_sio2,           0.0,                          0.0   },  // a1
      {a_sio2 * (-0.5),  a_sio2 * 0.86602540378443865, 0.0   },  // a2
      {0.0,              0.0,                          c_sio2}};  // a3
  int npbc = 3;

  c_gfnff_calculator calc_pbc =
      c_gfnff_calculator_init_pbc(nat_pbc, at_pbc, xyz_pbc, 0, 1,
                                  lattice_sio2, npbc);

  if (calc_pbc.ptr == NULL) {
    printf("Error initializing PBC gfnff calculator.\n");
    return 1;
  }
  printf("PBC gfnff calculator initialized successfully.\n");

  double energy_pbc;
  double gradient_pbc[9][3];
  double sigma_pbc[3][3];  // stress tensor
  int iostat_pbc;

  c_gfnff_calculator_singlepoint(&calc_pbc, nat_pbc, at_pbc, xyz_pbc,
                                 &energy_pbc, gradient_pbc, sigma_pbc,
                                 lattice_sio2, &iostat_pbc);

  if (iostat_pbc == 0) {
    printf("PBC singlepoint calculation successful.\n");
    printf("PBC Energy: %f\n", energy_pbc);
    for (int i = 0; i < 3; i++) {
      int j = 0;
      printf("PBC Gradient[%d][%d] = %e\n", j, i, gradient_pbc[j][i]);
    }

    // Print the PBC stress tensor
    printf("Sigma (PBC):\n");
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        printf("  sigma_pbc[%d][%d] = %e\n", i, j, sigma_pbc[i][j]);
  } else {
    printf("PBC singlepoint calculation failed with iostat = %d\n", iostat_pbc);
  }

  // The Hessian must refuse a periodic calculator rather than quietly
  // returning a molecular one.
  {
    double hdummy[1] = {-12345.0};
    int h_iostat_pbc;
    c_gfnff_calculator_hessian(&calc_pbc, nat_pbc, at_pbc, xyz_pbc, hdummy,
                               NULL, NULL, 0.0, &h_iostat_pbc);
    if (h_iostat_pbc != GFNFF_C_PBC_UNSUPPORTED) {
      printf("PBC Hessian should have been rejected, got iostat = %d\n",
             h_iostat_pbc);
      return 1;
    }
    if (hdummy[0] != -12345.0) {
      printf("PBC Hessian wrote to the caller buffer despite failing.\n");
      return 1;
    }
    printf("PBC Hessian correctly rejected (iostat = %d, buffer untouched).\n",
           h_iostat_pbc);
  }

  c_gfnff_calculator_deallocate(&calc_pbc);
  if (calc_pbc.ptr == NULL) {
    printf("PBC gfnff calculator deallocated successfully.\n");
  }

  return 0;
}
