/** 
 * Here, we will use the precursor film model to simulate spreading of a drop on a superhydrophilic surface.
 * Microscopic contact sngle is 0 by definition.
 * Rest of the dynamics near the contact line is set to evolve based on the Navier-Stokes equations without any other assumption. 
 * We call this method: Precursor film + NS
*/

// 1 is drop and 2 is surrounding fluid
#include "axi.h"
#include "navier-stokes/centered.h"
#define FILTERED
#include "two-phase.h"
#include "tension.h"
#include "distance.h"
#include "adapt_wavelet_limited_v2.h"

#define MINlevel 3                                              // maximum level

#define tsnap (5e-4)

// Error tolerances
#define fErr (1e-3)                                 // error tolerance in VOF
#define KErr (1e-4)                                 // error tolerance in KAPPA
#define VelErr (1e-2)                            // error tolerances in velocity
#define OmegaErr (1e-2)                            // error tolerances in velocity

// boundary conditions
f[left] = dirichlet(1.0);
u.t[left] = dirichlet(0.0);


double Ohw, Oha, hf, tmax, Ldomain, delta;
int MAXlevel;

int main(int argc, char const *argv[]) {
  Ldomain = 2.1; // size of the domain
  delta = 1e-2; // initial bridge length connecting the drop to the substrate
  Ohw = 1e-2; // dimensionless viscosity of the spreading drop
  tmax = 1.0; // maximum simulation time
  MAXlevel = 11; // maximum level of refinement
  Oha = 1e-5; // dimensionless viscosity of the surrounding fluid
  hf = 1e-2; // height of the precursor film

  L0=Ldomain;
  X0=-hf*1.001; Y0=0.;
  init_grid (1 << (4));

  char comm[80];
  sprintf (comm, "mkdir -p intermediate");
  system(comm);

  rho1 = 1e0; mu1 = Ohw;
  rho2 = 1e-3; mu2 = Oha;

  f.sigma = 1.0;

  fprintf(ferr, "Level %d tmax %g. Ohw %3.2e, Oha %3.2e, hf %3.2f\n", MAXlevel, tmax, Ohw, Oha, hf);
  run();
}

int refRegion(double x, double y, double z){

  return (x < 1.0*hf && y < 2.0 ? MAXlevel:
          x < 10.0*hf && y < 2.5 ? MAXlevel-1:
          x < 1.0 && y < 3.0 ? MAXlevel-2:
          x < 2.0 && y < 4.0 ? MAXlevel-3:
          x < 4.0 && y < 4.0 ? MAXlevel-4:
          MAXlevel-5
        );
}

event init(t = 0){
  if(!restore (file = "dump")){
    char filename[60];
    /**
    Initialization for f
    */
    sprintf(filename, "f_init-%3.2f.dat", delta);

    FILE * fp = fopen(filename,"rb");
    if (fp == NULL){
      fprintf(ferr, "There is no file named %s\n", filename);
      return 1;
    }

    coord* InitialShape;
    InitialShape = input_xy(fp);
    fclose (fp);
    scalar d[];
    distance (d, InitialShape);
    while (adapt_wavelet_limited ((scalar *){f, d}, (double[]){1e-8, 1e-8}, refRegion).nf);
    /**
    The distance function is defined at the center of each cell, we have
    to calculate the value of this function at each vertex. */
    vertex scalar phi[];
    foreach_vertex(){
      phi[] = -(d[] + d[-1] + d[0,-1] + d[-1,-1])/4.;
    }
    /**
    We can now initialize the volume fractions in the domain. */
    fractions (phi, f);
  }
  dump (file = "dump");
  // return 1;
}

scalar KAPPA[], omega[];
event adapt(i++) {
  vorticity (u, omega);
  curvature(f, KAPPA);
  foreach(){
    omega[] *= f[];
  }
  adapt_wavelet_limited ((scalar *){f, u.x, u.y, KAPPA, omega},
    (double[]){fErr, VelErr, VelErr, KErr, OmegaErr},
    refRegion, MINlevel);
}

// Outputs
event writingFiles (t = 0; t += tsnap; t <= tmax + tsnap) {
  dump (file = "dump");
  char nameOut[80];
  sprintf (nameOut, "intermediate/snapshot-%5.4f", t);
  dump (file = nameOut);
}

event logWriting (i++) {
  double ke = 0.;
  foreach (reduction(+:ke)){
    ke += sq(Delta)*(sq(u.x[]) + sq(u.y[]))*rho(f[]);
  }
  static FILE * fp;
  if (i == 0) {
    fprintf (ferr, "i dt t ke\n");
    fp = fopen ("log", "w");
    fprintf (fp, "i dt t ke\n");
    fprintf (fp, "%d %g %g %g\n", i, dt, t, ke);
    fclose(fp);
  } else {
    fp = fopen ("log", "a");
    fprintf (fp, "%d %g %g %g\n", i, dt, t, ke);
    fclose(fp);
  }
  fprintf (ferr, "%d %g %g %g\n", i, dt, t, ke);
}
