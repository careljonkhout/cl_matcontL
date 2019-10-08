#include <nvector/nvector_serial.h>
#include <sundials/sundials_types.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunmatrix/sunmatrix_band.h>
#include <math.h>
#include <mex.h>
#include "user_data.h"


#define X_INDEX(i) (2*(i))
#define Y_INDEX(i) (2*(i)+1)

#define X(i) u[X_INDEX(i)]
#define Y(i) u[Y_INDEX(i)]

#if JACOBIAN_STORAGE == DENSE
#define JAC_XX(i,j) SM_ELEMENT_D(jacobian, X_INDEX(i), X_INDEX(j))
#define JAC_XY(i,j) SM_ELEMENT_D(jacobian, X_INDEX(i), Y_INDEX(j))
#define JAC_YX(i,j) SM_ELEMENT_D(jacobian, Y_INDEX(i), X_INDEX(j))
#define JAC_YY(i,j) SM_ELEMENT_D(jacobian, Y_INDEX(i), Y_INDEX(j))
#endif

#if JACOBIAN_STORAGE == BANDED
#define JAC_XX(i,j) SM_ELEMENT_B(jacobian, X_INDEX(i), X_INDEX(j))
#define JAC_XY(i,j) SM_ELEMENT_B(jacobian, X_INDEX(i), Y_INDEX(j))
#define JAC_YX(i,j) SM_ELEMENT_B(jacobian, Y_INDEX(i), X_INDEX(j))
#define JAC_YY(i,j) SM_ELEMENT_B(jacobian, Y_INDEX(i), Y_INDEX(j))
#endif

#define L  parameters[0]
#define A  parameters[1]
#define B  parameters[2]
#define DX parameters[3]
#define DY parameters[4]

int jacobian_dydt(
               realtype t, N_Vector u_vec, N_Vector fy, SUNMatrix jacobian,
              void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
  
  realtype* u   = N_VGetArrayPointer(u_vec);
  realtype* parameters = ((UserData) user_data)->parameters;
  
  double cx = DX * (N_MESH_POINTS + 1) * (N_MESH_POINTS+1) / (L*L);
  double cy = DY * (N_MESH_POINTS + 1) * (N_MESH_POINTS+1) / (L*L);
  
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
    JAC_XX(i,i) = - 2 * cx - (B+1) + 2 * X(i) * Y(i);
    JAC_XY(i,i) =   X(i) * X(i);
    JAC_YX(i,i) =   B - 2 * X(i) * Y(i);
    JAC_YY(i,i) = - 2 * cy - X(i) * X(i);
  }
  
  for ( int i = 0; i < N_MESH_POINTS - 1; i++ ) {
    JAC_XX(i    , i + 1) = cx;
    JAC_XX(i + 1, i    ) = cx;
    JAC_YY(i    , i + 1) = cy;
    JAC_YY(i + 1, i    ) = cy;
  }
  return 0;
}