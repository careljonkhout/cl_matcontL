// defines the Jacobian matrix of the Elezgaray-Aneodo system

#include <nvector/nvector_serial.h>
#include <sundials/sundials_types.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunmatrix/sunmatrix_band.h>
#include <math.h>
#include <mex.h>
#include "user_data.h"


#define U_INDEX(i) (2*(i))
#define V_INDEX(i) (2*(i)+1)

#define U(i) y[U_INDEX(i)]
#define V(i) y[V_INDEX(i)]

#if JACOBIAN_STORAGE == DENSE
#define JAC_UU(i,j) SM_ELEMENT_D(jacobian, U_INDEX(i), U_INDEX(j))
#define JAC_UV(i,j) SM_ELEMENT_D(jacobian, U_INDEX(i), V_INDEX(j))
#define JAC_VU(i,j) SM_ELEMENT_D(jacobian, V_INDEX(i), U_INDEX(j))
#define JAC_VV(i,j) SM_ELEMENT_D(jacobian, V_INDEX(i), V_INDEX(j))
#endif

#if JACOBIAN_STORAGE == BANDED
#define JAC_UU(i,j) SM_ELEMENT_B(jacobian, U_INDEX(i), U_INDEX(j))
#define JAC_UV(i,j) SM_ELEMENT_B(jacobian, U_INDEX(i), V_INDEX(j))
#define JAC_VU(i,j) SM_ELEMENT_B(jacobian, V_INDEX(i), U_INDEX(j))
#define JAC_VV(i,j) SM_ELEMENT_B(jacobian, V_INDEX(i), V_INDEX(j))
#endif

#define D     parameters[0]
#define eps   parameters[1]
#define alpha parameters[2]

int jacobian_dydt(
               realtype t, N_Vector y_vec, N_Vector fy, SUNMatrix jacobian,
              void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
  
  realtype* y   = N_VGetArrayPointer(y_vec);
  realtype* parameters = ((UserData) user_data)->parameters;
  
  double c = D * (N_MESH_POINTS + 1) * (N_MESH_POINTS + 1);
  
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
  for ( int i = 0; i < N_MESH_POINTS; i++ ) {
    JAC_UU(i,i) = - 2 * c - (2 * U(i)  + 3 * U(i) * U(i)) / eps;
    JAC_UV(i,i) =   1 / eps;
    JAC_VU(i,i) = - 1;
    JAC_VV(i,i) = - 2 * c;
  }
  
  for ( int i = 0; i < N_MESH_POINTS - 1; i++ ) {
    JAC_UU(i    , i + 1) = c;
    JAC_UU(i + 1, i    ) = c;
    JAC_VV(i    , i + 1) = c;
    JAC_VV(i + 1, i    ) = c;
  }
  return 0;
}