#include <nvector/nvector_serial.h>
#include <sundials/sundials_types.h>
#include <math.h>
#include "user_data.h"

int d_sensitivity_dt_pars(int Ns, realtype t,
                N_Vector y_vector, 
                N_Vector ydot,
                int parameter_index,
                N_Vector s_vector, 
                N_Vector ds_vector,
                void *user_data,
                N_Vector tmp1, N_Vector tmp2) {
   
  realtype* y = N_VGetArrayPointer(y_vector);
                
    
  realtype* s          = N_VGetArrayPointer( s_vector);
  realtype* ds         = N_VGetArrayPointer(ds_vector);
  realtype* parameters = ((UserData) user_data)->parameters;
 
  switch (parameter_index) {
  case 0:
  ds[0] =  -parameters[2]*s[0]*exp(-y[0])-parameters[1]*parameters[3]*s[3]*exp(y[2])-parameters[1]*s[2]*exp(y[2])*(parameters[3]*y[3]+1.0);
  ds[1] =  parameters[1]*parameters[3]*s[3]*exp(-y[1]+y[2]+y[0])+parameters[1]*s[0]*exp(-y[1]+y[2]+y[0])*(parameters[3]*y[3]+1.0)-parameters[1]*s[1]*exp(-y[1]+y[2]+y[0])*(parameters[3]*y[3]+1.0)+parameters[1]*s[2]*exp(-y[1]+y[2]+y[0])*(parameters[3]*y[3]+1.0)-1.0;
  ds[2] =  exp(y[1]-y[2])+parameters[0]*s[1]*exp(y[1]-y[2])-parameters[0]*s[2]*exp(y[1]-y[2]);
  ds[3] =  -s[4]*(3.141592653589793*2.0+y[3]*y[4]*2.0)-s[3]*((y[3]*y[3])*3.0+y[4]*y[4]-1.0);
  ds[4] =  s[3]*(3.141592653589793*2.0-y[3]*y[4]*2.0)-s[4]*(y[3]*y[3]+(y[4]*y[4])*3.0-1.0);
  break;
  case 1:
  ds[0] =  -exp(y[2])*(parameters[3]*y[3]+1.0)-parameters[2]*s[0]*exp(-y[0])-parameters[1]*parameters[3]*s[3]*exp(y[2])-parameters[1]*s[2]*exp(y[2])*(parameters[3]*y[3]+1.0);
  ds[1] =  exp(-y[1]+y[2]+y[0])*(parameters[3]*y[3]+1.0)+parameters[1]*parameters[3]*s[3]*exp(-y[1]+y[2]+y[0])+parameters[1]*s[0]*exp(-y[1]+y[2]+y[0])*(parameters[3]*y[3]+1.0)-parameters[1]*s[1]*exp(-y[1]+y[2]+y[0])*(parameters[3]*y[3]+1.0)+parameters[1]*s[2]*exp(-y[1]+y[2]+y[0])*(parameters[3]*y[3]+1.0);
  ds[2] =  parameters[0]*s[1]*exp(y[1]-y[2])-parameters[0]*s[2]*exp(y[1]-y[2]);
  ds[3] =  -s[4]*(3.141592653589793*2.0+y[3]*y[4]*2.0)-s[3]*((y[3]*y[3])*3.0+y[4]*y[4]-1.0);
  ds[4] =  s[3]*(3.141592653589793*2.0-y[3]*y[4]*2.0)-s[4]*(y[3]*y[3]+(y[4]*y[4])*3.0-1.0);
  break;
  case 2:
  ds[0] =  exp(-y[0])-parameters[2]*s[0]*exp(-y[0])-parameters[1]*parameters[3]*s[3]*exp(y[2])-parameters[1]*s[2]*exp(y[2])*(parameters[3]*y[3]+1.0)-1.0;
  ds[1] =  parameters[1]*parameters[3]*s[3]*exp(-y[1]+y[2]+y[0])+parameters[1]*s[0]*exp(-y[1]+y[2]+y[0])*(parameters[3]*y[3]+1.0)-parameters[1]*s[1]*exp(-y[1]+y[2]+y[0])*(parameters[3]*y[3]+1.0)+parameters[1]*s[2]*exp(-y[1]+y[2]+y[0])*(parameters[3]*y[3]+1.0)-1.0;
  ds[2] =  parameters[0]*s[1]*exp(y[1]-y[2])-parameters[0]*s[2]*exp(y[1]-y[2])-1.0;
  ds[3] =  -s[4]*(3.141592653589793*2.0+y[3]*y[4]*2.0)-s[3]*((y[3]*y[3])*3.0+y[4]*y[4]-1.0);
  ds[4] =  s[3]*(3.141592653589793*2.0-y[3]*y[4]*2.0)-s[4]*(y[3]*y[3]+(y[4]*y[4])*3.0-1.0);
  break;
  case 3:
  ds[0] =  -parameters[1]*y[3]*exp(y[2])-parameters[2]*s[0]*exp(-y[0])-parameters[1]*parameters[3]*s[3]*exp(y[2])-parameters[1]*s[2]*exp(y[2])*(parameters[3]*y[3]+1.0);
  ds[1] =  parameters[1]*y[3]*exp(-y[1]+y[2]+y[0])+parameters[1]*parameters[3]*s[3]*exp(-y[1]+y[2]+y[0])+parameters[1]*s[0]*exp(-y[1]+y[2]+y[0])*(parameters[3]*y[3]+1.0)-parameters[1]*s[1]*exp(-y[1]+y[2]+y[0])*(parameters[3]*y[3]+1.0)+parameters[1]*s[2]*exp(-y[1]+y[2]+y[0])*(parameters[3]*y[3]+1.0);
  ds[2] =  parameters[0]*s[1]*exp(y[1]-y[2])-parameters[0]*s[2]*exp(y[1]-y[2]);
  ds[3] =  -s[4]*(3.141592653589793*2.0+y[3]*y[4]*2.0)-s[3]*((y[3]*y[3])*3.0+y[4]*y[4]-1.0);
  ds[4] =  s[3]*(3.141592653589793*2.0-y[3]*y[4]*2.0)-s[4]*(y[3]*y[3]+(y[4]*y[4])*3.0-1.0);
  break;
  case 4:
  ds[0] =  -parameters[2]*s[0]*exp(-y[0])-parameters[1]*parameters[3]*s[3]*exp(y[2])-parameters[1]*s[2]*exp(y[2])*(parameters[3]*y[3]+1.0);
  ds[1] =  parameters[1]*parameters[3]*s[3]*exp(-y[1]+y[2]+y[0])+parameters[1]*s[0]*exp(-y[1]+y[2]+y[0])*(parameters[3]*y[3]+1.0)-parameters[1]*s[1]*exp(-y[1]+y[2]+y[0])*(parameters[3]*y[3]+1.0)+parameters[1]*s[2]*exp(-y[1]+y[2]+y[0])*(parameters[3]*y[3]+1.0);
  ds[2] =  parameters[0]*s[1]*exp(y[1]-y[2])-parameters[0]*s[2]*exp(y[1]-y[2])-1.0;
  ds[3] =  -s[4]*(3.141592653589793*2.0+y[3]*y[4]*2.0)-s[3]*((y[3]*y[3])*3.0+y[4]*y[4]-1.0);
  ds[4] =  s[3]*(3.141592653589793*2.0-y[3]*y[4]*2.0)-s[4]*(y[3]*y[3]+(y[4]*y[4])*3.0-1.0);
  break;
  }


  return 0;
}