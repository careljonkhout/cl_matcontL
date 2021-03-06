// this file was generated by the system file generator of cl_matcontL

// this mex function computes the Jacobian matrix of system of ODEs called
//                       "brusselator_1d_N_50"

#include <math.h>
#include "mex.h"
#include "user_data.h"


#define N_INPUTS     N_PARAMETERS + 2

#define INPUT_T 0
#define INPUT_Y 1


void jacobian_dydt_mex(double* u, double* parameters, double* jacobian);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

  double parameters[N_PARAMETERS];  
  double *y;                      /* 1xN input matrix */
  double *outMatrix;              /* output matrix */

  /* check for proper number of arguments */
  if(nrhs != N_INPUTS) {
    mexErrMsgIdAndTxt("brusselator_1d_N_50_dydt:nrhs","Five inputs required.");
  }

  if ( !mxIsDouble(prhs[INPUT_Y]) ) {
    mexErrMsgIdAndTxt("system_brusselator_1d_N_50:not_double",
                      "Error: Input vector y is not a double.");
  }

  if ( mxGetNumberOfElements(prhs[INPUT_Y]) != NEQ ) {
    mexErrMsgIdAndTxt("system_brusselator_1d_N_50:wrong_size",
                     "Input vector y must have 100 elements.");
  }

  for ( int i = 0; i < N_PARAMETERS; i++ ) {
    if ( !mxIsDouble(prhs[i + 2]) ) {
      mexErrMsgIdAndTxt("brusselator_1d_N_50_dydt:not_double",
                        "Error: one of the parameters is not a double");
    } else if ( mxIsComplex(prhs[i + 2]) ) {
      mexErrMsgIdAndTxt("brusselator_1d_N_50_dydt:not_complex",
                        "Error: one of the parameters not a real number");
    } else if ( mxGetNumberOfElements(prhs[i + 2]) != 1 ) {
      mexErrMsgIdAndTxt("brusselator_1d_N_50_dydt:notScalar",
                        "Error: one of the parameters is not scalar");
    } else {
      parameters[i] = mxGetScalar(prhs[i+2]);
    }
  }
    
    
  /* create a pointer to the real data in the input matrix  */
  #if MX_HAS_INTERLEAVED_COMPLEX
  y = mxGetDoubles(prhs[INPUT_Y]);
  #else
  y = mxGetPr(prhs[INPUT_Y]);
  #endif

  /* get dimensions of the input matrix */

  /* create the output matrix */
  plhs[0] = mxCreateDoubleMatrix(NEQ, NEQ, mxREAL);

  /* get a pointer to the real data in the output matrix */
  #if MX_HAS_INTERLEAVED_COMPLEX
  outMatrix = mxGetDoubles(plhs[0]);
  #else
  outMatrix = mxGetPr(plhs[0]);
  #endif

  /* call the computational routine */
  jacobian_dydt_mex(y,parameters,outMatrix);
    
}

#define U_INDEX(i) (2*(i))
#define V_INDEX(i) (2*(i)+1)

#define U(i) y[U_INDEX(i)]
#define V(i) y[V_INDEX(i)]


#define ELEMENT(mat,i,j) mat[i + j*NEQ]

#define JAC_UU(i,j) ELEMENT(jacobian, U_INDEX(i), U_INDEX(j))
#define JAC_UV(i,j) ELEMENT(jacobian, U_INDEX(i), V_INDEX(j))
#define JAC_VU(i,j) ELEMENT(jacobian, V_INDEX(i), U_INDEX(j))
#define JAC_VV(i,j) ELEMENT(jacobian, V_INDEX(i), V_INDEX(j))

#define D     parameters[0]
#define eps   parameters[1]
#define alpha parameters[2]


void jacobian_dydt_mex(double* y, double* parameters, double* jacobian) {
  
  double c = D * (N_MESH_POINTS + 1) * (N_MESH_POINTS + 1);
  
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
}