#include <nvector/nvector_serial.h>
#include <sundials/sundials_types.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <math.h>
#include "user_data.h"

int jacobian_dydt(
               realtype t, N_Vector y_vec, N_Vector fy, SUNMatrix jac_structure,
              void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
  
  realtype* y   = N_VGetArrayPointer(y_vec);
  realtype* jac = SM_DATA_D(jac_structure); // get the pointer to the jac data
  realtype* parameters = ((UserData) user_data)->parameters;

  jac[0] =  parameters[0];
  jac[1] =  -parameters[1];
  jac[2] =  y[2];
  jac[3] =  -y[1];
  jac[4] =  -parameters[0];
  jac[6] =  1.0;
  jac[7] =  -y[0];
  jac[10] =  y[0];
  jac[11] =  parameters[2];
  
  return 0;
}