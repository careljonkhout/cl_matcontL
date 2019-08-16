// this file was generated by the system file generator of cl_matcontL

// this mex function computes the Jacobian matrix of system of ODEs called
//                       "lorenz"

#include <math.h>
#include "mex.h"

#define n_parameters 3
#define n_inputs     5
#define n_outputs    1
#define system_size  3

#define input_t 0
#define input_y 1


void lorenz_jacobian(double* y, double* parameters, double* dydt);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

  double parameters[n_parameters];  
  double *y;                      /* 1xN input matrix */
  size_t ncols;                   /* size of matrix */
  double *outMatrix;              /* output matrix */

  /* check for proper number of arguments */
  if(nrhs != n_inputs) {
    mexErrMsgIdAndTxt("lorenz_dydt:nrhs","Five inputs required.");
  }
  if(nlhs != n_outputs) {
    mexErrMsgIdAndTxt("lorenz_dydt:nlhs","One output required.");
  }


  if ( !mxIsDouble(prhs[input_y]) ) {
    mexErrMsgIdAndTxt("system_lorenz:not_double",
                      "Error: Input vector y is not a double.");
  }

  if ( mxGetNumberOfElements(prhs[input_y]) != system_size ) {
    mexErrMsgIdAndTxt("system_lorenz:wrong_size",
                     "Input vector y must have 3 elements.");
  }

  for ( int i = 0; i < n_parameters; i++ ) {
    if ( !mxIsDouble(prhs[i + 2]) ) {
      mexErrMsgIdAndTxt("lorenz_dydt:not_double",
                        "Error: one of the parameters is not a double");
    } else if ( mxIsComplex(prhs[i + 2]) ) {
      mexErrMsgIdAndTxt("lorenz_dydt:not_complex",
                        "Error: one of the parameters not a real number");
    } else if ( mxGetNumberOfElements(prhs[i + 2]) != 1 ) {
      mexErrMsgIdAndTxt("lorenz_dydt:notScalar",
                        "Error: one of the parameters is not scalar");
    } else {
      parameters[i] = mxGetScalar(prhs[i+2]);
    }
  }
    
    
  /* create a pointer to the real data in the input matrix  */
  #if MX_HAS_INTERLEAVED_COMPLEX
  y = mxGetDoubles(prhs[input_y]);
  #else
  y = mxGetPr(prhs[input_y]);
  #endif

  /* get dimensions of the input matrix */

  /* create the output matrix */
  plhs[0] = mxCreateDoubleMatrix(system_size, system_size, mxREAL);

  /* get a pointer to the real data in the output matrix */
  #if MX_HAS_INTERLEAVED_COMPLEX
  outMatrix = mxGetDoubles(plhs[0]);
  #else
  outMatrix = mxGetPr(plhs[0]);
  #endif

  /* call the computational routine */
  lorenz_jacobian(y,parameters,outMatrix);
    
}

void lorenz_jacobian(double* y, double* parameters, double* jac) {
  jac[0] =  -parameters[0];
  jac[1] =  parameters[1]-y[2];
  jac[2] =  y[1];
  jac[3] =  parameters[0];
  jac[4] =  -1.0;
  jac[5] =  y[0];
  jac[7] =  -y[0];
  jac[8] =  -parameters[2];
}