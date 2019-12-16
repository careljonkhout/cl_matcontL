#include <nvector/nvector_serial.h>
#include <sundials/sundials_types.h>
#include <math.h>
#include <mex.h>
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
                
  mexErrMsgIdAndTxt("cvode:not_implemented", 
                    "d_sensitivity_dt_pars is not implemented");
  return 0;
}