// defines the Elezgaray-Aneodo system

#include <math.h>
#include "user_data.h"
#include <mex.h>

#define u(i) y[2*(i)]
#define v(i) y[2*(i)+1]

#define dudt(i) dydt[2*(i)]
#define dvdt(i) dydt[2*(i)+1]

#define D     parameters[0]
#define eps   parameters[1]
#define alpha parameters[2]

#define reaction_u(i) (v(i) - u(i) * u(i) * (u(i) + 1)) / eps
#define reaction_v(i)  alpha - u(i)


void dydt_ode(double* y, double* dydt, double* parameters) {
  
  // the next four lines define the values of the Dirichlet boundary conditions
  double u0 = -2;
  double u1 = -2;
  double v0 = -4;
  double v1 = -4;
   
  double c = D * (N_MESH_POINTS + 1) * (N_MESH_POINTS + 1);
  
  dudt(0) = ( u0 - 2 * u(0) + u(1) ) * c + reaction_u(0);
  dvdt(0) = ( v0 - 2 * v(0) + v(1) ) * c + reaction_v(0);
  
  // #pragma omp parallel for
  for (int i = 1; i < N_MESH_POINTS - 1; i++) {
    dudt(i) = ( u(i-1) - 2 * u(i) + u(i+1) ) * c + reaction_u(i);
    dvdt(i) = ( v(i-1) - 2 * v(i) + v(i+1) ) * c + reaction_v(i);
  }
  
  int i = N_MESH_POINTS - 1;
  
  dudt(i) = ( u(i-1) - 2 * u(i) + u1 ) * c + reaction_u(i);
  dvdt(i) = ( v(i-1) - 2 * v(i) + v1 ) * c + reaction_v(i);
  
}
