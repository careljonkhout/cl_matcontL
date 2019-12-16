#include <nvector/nvector_serial.h>
#include <sundials/sundials_types.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunmatrix/sunmatrix_band.h>
#include <math.h>
#include <mex.h>
#include "user_data.h"

#if JACOBIAN_STORAGE == DENSE
#define JAC(i,j) (SM_ELEMENT_D(jacobian, (i), (j)))
#endif

#if JACOBIAN_STORAGE == BANDED
#define JAC(i,j) (SM_ELEMENT_B(jacobian, (i), (j)))
#endif

int jacobian_dydt(
               realtype t, N_Vector y_vec, N_Vector fy, SUNMatrix jacobian,
              void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
  
  realtype* y          = N_VGetArrayPointer(y_vec);
  realtype* parameters = ((UserData) user_data)->parameters;

  #if JACOBIAN_STORAGE == DENSE
  if (SUNMatGetID(jacobian) != SUNMATRIX_DENSE) {
    mexErrMsgIdAndTxt("cvode:wrong_matrix_type",
            "expected type %d = SUNMATRIX_DENSE, got type %d",
            SUNMATRIX_DENSE, SUNMatGetID(jacobian));
  }
  #endif
  
  #if JACOBIAN_STORAGE == BANDED
  if (SUNMatGetID(jacobian) != SUNMATRIX_BAND) {
    mexErrMsgIdAndTxt("cvode:wrong_matrix_type",
            "expected type %d = SUNMATRIX_BAND, got type %d",
            SUNMATRIX_DENSE, SUNMatGetID(jacobian));
  }
  #endif

  JAC(0,0) =  -parameters[2]*exp(-y[0]);
  JAC(1,0) =  parameters[1]*exp(-y[1]+y[2]+y[0])*(parameters[3]*y[3]+1.0);
  JAC(1,1) =  -parameters[1]*exp(-y[1]+y[2]+y[0])*(parameters[3]*y[3]+1.0);
  JAC(2,1) =  parameters[0]*exp(y[1]-y[2]);
  JAC(0,2) =  -parameters[1]*exp(y[2])*(parameters[3]*y[3]+1.0);
  JAC(1,2) =  parameters[1]*exp(-y[1]+y[2]+y[0])*(parameters[3]*y[3]+1.0);
  JAC(2,2) =  -parameters[0]*exp(y[1]-y[2]);
  JAC(0,3) =  -parameters[1]*parameters[3]*exp(y[2]);
  JAC(1,3) =  parameters[1]*parameters[3]*exp(-y[1]+y[2]+y[0]);
  JAC(3,3) =  (y[3]*y[3])*-3.0-y[4]*y[4]+1.0;
  JAC(4,3) =  3.141592653589793*2.0-y[3]*y[4]*2.0;
  JAC(3,4) =  3.141592653589793*-2.0-y[3]*y[4]*2.0;
  JAC(4,4) =  -y[3]*y[3]-(y[4]*y[4])*3.0+1.0;
  
  return 0;
}