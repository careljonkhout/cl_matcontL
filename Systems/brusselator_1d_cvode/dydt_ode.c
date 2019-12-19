#include <math.h>
#include "user_data.h"
#include <mex.h>

#define x(i) u[2*(i)]
#define y(i) u[2*(i)+1]

#define dxdt(i) dudt[2*(i)]
#define dydt(i) dudt[2*(i)+1]

#define L  parameters[0]
#define A  parameters[1]
#define B  parameters[2]
#define DX parameters[3]
#define DY parameters[4]

#define reaction_x(i) A   -   (B + 1) * x(i)   +   x(i) * x(i) * y(i)
#define reaction_y(i)          B      * x(i)   -   x(i) * x(i) * y(i)


void dydt_ode(double* u, double* dudt, double* parameters) {
  
  // the next four lines define the values of the Dirichlet boundary conditions
  double x0 = A;
  double x1 = A;
  double y0 = B/A;
  double y1 = B/A;
   
  double cx = DX * (N_MESH_POINTS + 1) * (N_MESH_POINTS + 1) / (L * L);
  double cy = DY * (N_MESH_POINTS + 1) * (N_MESH_POINTS + 1) / (L * L);
  
  dxdt(0) = ( x0 - 2 * x(0) + x(1) ) * cx + reaction_x(0);
  dydt(0) = ( y0 - 2 * y(0) + y(1) ) * cy + reaction_y(0);
  
  for (int i = 1; i < N_MESH_POINTS - 1; i++) {
    dxdt(i) = ( x(i-1) - 2 * x(i) + x(i+1) ) * cx + reaction_x(i);
    dydt(i) = ( y(i-1) - 2 * y(i) + y(i+1) ) * cy + reaction_y(i);
  }
  
  int i = N_MESH_POINTS - 1;
  
  dxdt(i) = ( x(i-1) - 2 * x(i) + x1 ) * cx + reaction_x(i);
  dydt(i) = ( y(i-1) - 2 * y(i) + y1 ) * cy + reaction_y(i);
  
}
