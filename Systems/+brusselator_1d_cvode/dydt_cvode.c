#include <nvector/nvector_serial.h>
#include <sundials/sundials_types.h>
#include <math.h>
#include "user_data.h"

void dydt_ode(double* y, double* dydt, double* parameters);

int dydt_cvode(realtype t, N_Vector y_vec, N_Vector dydt_vec, void *user_data) {
  
  realtype* y              = N_VGetArrayPointer(y_vec);
  realtype* dydt           = N_VGetArrayPointer(dydt_vec);
  realtype* parameters     = ((UserData) user_data)->parameters;
  
  dydt_ode(y,dydt,parameters);
  
  return 0;
}
