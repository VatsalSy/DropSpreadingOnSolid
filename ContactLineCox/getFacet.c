/* Title: Getting Facets
# Author: Vatsal Sanjay
# vatsalsanjay@gmail.com
# Physics of Fluids
*/
#include "navier-stokes/centered.h"
scalar f[], * interfaces = {f};
#include "contact.h"
vector h[];
double theta0, M_c;

char filename[80];
int main(int a, char const *arguments[])
{
  theta0 = 15.0;
  M_c = atof(arguments[2]);
  // boundary conditions
  h.t[left] = contact_angle (pow((pow(theta0*pi/180., 3) + M_c*u.t[]), 1./3.));
  u.t[left] = dirichlet(0.0);
  
  sprintf (filename, "%s", arguments[1]);
  restore (file = filename);
  #if TREE
    f.prolongation = fraction_refine;
  #endif
  boundary((scalar *){f});

  FILE * fp = ferr;
  output_facets(f,fp);
  fflush (fp);
  fclose (fp);
}
