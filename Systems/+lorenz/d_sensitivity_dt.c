#include <nvector/nvector_serial.h>
#include <sundials/sundials_types.h>
#include <math.h>
#include "user_data.h"

int d_sensitivity_dt(int Ns, realtype t,
                N_Vector y_vector, 
                N_Vector ydot,
                int iS,
                N_Vector s_vector, 
                N_Vector ds_vector,
                void *user_data,
                N_Vector tmp1, N_Vector tmp2) {
   
  double* y = N_VGetArrayPointer(y_vector);
                
    
  double* s  = N_VGetArrayPointer( s_vector);
  double* ds = N_VGetArrayPointer(ds_vector);
  double* parameters = ((UserData) user_data)->parameters;
 
  ds[0] =  parameters[0]*s[0]-parameters[0]*s[1];
  ds[1] =  -parameters[1]*s[0];
  ds[2] =  s[1]+s[2]*y[0]+s[0]*y[2];
  ds[3] =  parameters[2]*s[2]-s[1]*y[0]-s[0]*y[1];

  return 0;
}