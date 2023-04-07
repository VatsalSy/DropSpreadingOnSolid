/* Tracking footprint of a spreading bubble
# Last Update: November 12, 2020
# Author: Vatsal Sanjay
# vatsalsanjay@gmail.com
# Physics of Fluids
*/
#include "navier-stokes/centered.h"
scalar f[], * interfaces = {f};
#include "contact.h"
vector h[];

char filename[80], nameTrack[80];
double DistCutoff;
double wt, xTP, yTP, vTP;
double theta0, M_c;

int main(int a, char const *arguments[])
{
  theta0 = 15.0;
  M_c = atof(arguments[3]);

  // boundary conditions
  h.t[left] = contact_angle (pow((pow(theta0*pi/180., 3) + M_c*u.t[]), 1./3.));
  u.t[left] = dirichlet(0.0);

  sprintf (filename, "%s", arguments[1]);
  sprintf(nameTrack, "%s", arguments[2]);

  DistCutoff = 1e-2;

  restore (file = filename);
  f.prolongation = fraction_refine;
  boundary((scalar *){f});

  // double x0 = -4., y0 = 0.;
  face vector s[];
  s.x.i = -1;
  wt = 0., xTP = 0., yTP = 0., vTP = 0.;
  foreach(){
    if (f[] > 1e-6 && f[] < 1. - 1e-6 && x > DistCutoff && x < 2*DistCutoff) {
      coord n = facet_normal (point, f, s);
      double alpha = plane_alpha (f[], n);
      coord segment[2];
      if (facets (n, alpha, segment) == 2){
        double x1 = x + (segment[0].x+segment[1].x)*Delta/2.;
        double y1 = y + (segment[0].y+segment[1].y)*Delta/2.;
        xTP += x1*Delta;
        yTP += y1*Delta;
        vTP += u.y[]*Delta;
        wt += Delta;
        }
      }
  }
  if (wt > 0.){
    xTP /= wt;
    yTP /= wt;
    vTP /= wt;
  } 

  fprintf(ferr, "%f %f %f %4.3e\n", t, xTP, yTP, vTP);
  FILE *fp2;
  fp2 = fopen (nameTrack, "a");
  fprintf(fp2, "%f %f %f %4.3e\n", t, xTP, yTP, vTP);
  fclose(fp2);
  // fprintf(ferr, "%f %f\n", x0, y0);
}
