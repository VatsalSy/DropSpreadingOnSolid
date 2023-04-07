/** 
 * Here, we simulate the spreading of a drop on a partially wetting surface.
 * Microscopic contact angle is set in the first grid cell at the left boundary.
 * Rest of the dynamics near the contact line is set to evolve based on the Navier-Stokes equations without any other assumption. 
 * We call this method: Contact line + NS
*/

// 1 is drop and 2 is surrounding fluid
#include "navier-stokes/centered.h"
#include "contact.h"
#define FILTERED
#include "two-phase.h"
#include "tension.h"
#include "distance.h"
#include "adapt_wavelet_limited_v2.h"
vector h[];

#define MINlevel 3                                              // maximum level

#define tsnap (5e-4)

// Error tolerances
#define fErr (1e-3)                                 // error tolerance in VOF
#define KErr (1e-4)                                 // error tolerance in KAPPA
#define VelErr (1e-2)                            // error tolerances in velocity
#define OmegaErr (1e-2)                            // error tolerances in velocity


double Ohw, Oha, hf, tmax, Ldomain, delta, theta0;
int MAXlevel;

h.t[left] = contact_angle(theta0*pi/180.);
u.t[left] = dirichlet(0.0);

int main(int argc, char const *argv[]) {
  
  Ldomain = 2.1; // size of the domain
  delta = 1e-2; // initial bridge length connecting the drop to the substrate
  Ohw = 1e-2; // dimensionless viscosity of the spreading drop
  tmax = 1.0; // maximum simulation time
  MAXlevel = 11; // maximum level of refinement
  Oha = 1e-5; // dimensionless viscosity of the surrounding fluid
  theta0 = 15.0; // contact angle of the substrate

  L0=Ldomain;
  X0=delta/8.; Y0=0.;
  init_grid (1 << (4));

  char comm[80];
  sprintf (comm, "mkdir -p intermediate");
  system(comm);

  rho1 = 1e0; mu1 = Ohw;
  rho2 = 1e-3; mu2 = Oha;

  f.height = h;
  f.sigma = 1.0;

  fprintf(ferr, "Level %d tmax %g. Ohw %3.2e, Oha %3.2e\n", MAXlevel, tmax, Ohw, Oha);
  run();
}

int refRegion(double x, double y, double z){

  return (x < 2.5*delta && y < 1.0 ? MAXlevel+1:
          x < 5.0*delta && y < 1.5 ? MAXlevel:
          x < 10*delta && y < 1.5 ? MAXlevel-1:
          x < 1.0 && y < 2.0 ? MAXlevel-2:
          x < 3.0 && y < 3.0 ? MAXlevel-3:
          MAXlevel-4
        );
}

event init(t = 0){

  if(!restore (file = "dump")){
    char filename1[60];
    /**
    Initialization for f
    */
    sprintf(filename1, "f1_init-%3.2f.dat", delta);

    FILE * fp1 = fopen(filename1,"rb");
    if (fp1 == NULL){
      fprintf(ferr, "There is no file named %s\n", filename1);
      return 1;
    }
    coord* InitialShape1;
    InitialShape1 = input_xy(fp1);
    fclose (fp1);
    scalar d1[];
    distance (d1, InitialShape1);
    while (adapt_wavelet_limited ((scalar *){f, d1}, (double[]){1e-8, 1e-8}, refRegion).nf);
    /**
    The distance function is defined at the center of each cell, we have
    to calculate the value of this function at each vertex. */
    vertex scalar phi1[];
    foreach_vertex(){
      phi1[] = -(d1[] + d1[-1] + d1[0,-1] + d1[-1,-1])/4.;
    }
    /**
    We can now initialize the volume fractions in the domain. */
    fractions (phi1, f);
  }
  dump (file = "dump");
  // return 1;
}

scalar KAPPA1[], omega[];
event adapt(i++) {
  vorticity (u, omega);
  curvature(f, KAPPA1);
  foreach(){
    omega[] *= f[];
  }
  adapt_wavelet_limited ((scalar *){f, u.x, u.y, KAPPA1, omega},
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
    ke += sq(Delta)*(sq(u.x[]) + sq(u.y[]))*(f[]);
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
