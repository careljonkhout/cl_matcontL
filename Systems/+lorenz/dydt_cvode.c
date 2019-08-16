#include <nvector/nvector_serial.h>
#include <sundials/sundials_types.h>
#include <math.h>
#include "user_data.h"

int dydt_ode(realtype t, N_Vector y_vec, N_Vector dydt_vec, void *user_data) {
  
  realtype* y    = N_VGetArrayPointer(y_vec);
  realtype* dydt = N_VGetArrayPointer(dydt_vec);
  realtype* parameters     = ((UserData) user_data)->parameters;

  dydt[0] =  -parameters[0]*(y[0]-y[1]);
  dydt[1] =  -y[1]+parameters[1]*y[0]-y[0]*y[2];
  dydt[2] =  -parameters[2]*y[2]+y[0]*y[1];

  return 0;
}
