#include <nvector/nvector_serial.h>
#include <sundials/sundials_types.h>
#include <math.h>
#include "user_data.h"

int dydt_cvode(realtype t, N_Vector y_vec, N_Vector dydt_vec, void *user_data) {
  
  realtype* y    = N_VGetArrayPointer(y_vec);
  realtype* dydt = N_VGetArrayPointer(dydt_vec);
  realtype* parameters     = ((UserData) user_data)->parameters;

  dydt[0] =  -parameters[2]+parameters[2]*exp(-y[0])-parameters[1]*exp(y[2])*(parameters[3]*y[3]+1.0);
  dydt[1] =  -parameters[0]-parameters[2]+parameters[1]*exp(-y[1]+y[2]+y[0])*(parameters[3]*y[3]+1.0);
  dydt[2] =  -parameters[4]-parameters[2]+parameters[0]*exp(y[1]-y[2]);
  dydt[3] =  y[3]-y[4]*3.141592653589793*2.0-y[3]*(y[3]*y[3]+y[4]*y[4]);
  dydt[4] =  y[4]+y[3]*3.141592653589793*2.0-y[4]*(y[3]*y[3]+y[4]*y[4]);

  return 0;
}
