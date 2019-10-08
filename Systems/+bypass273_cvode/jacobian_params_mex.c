// this file was generated by the system file generator of cl_matcontL

// this mex function computes the Jacobian matrix of system of ODEs called
//                       "bypass273_cvode"

#include <math.h>
#include "mex.h"

#define n_parameters 1
#define n_inputs     3
#define n_outputs    1
#define system_size  273

#define input_t 0
#define input_y 1
#define checks false

void bypass273_cvode_jacobian_params(double* y, double* parameters, double* dydt);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

  double parameters[n_parameters];  
  double *y;                      /* 1xN input matrix */
  size_t ncols;                   /* size of matrix */
  double *outMatrix;              /* output matrix */

  /* check for proper number of arguments */
  if(nrhs != n_inputs) {
    mexErrMsgIdAndTxt("bypass273_cvode_dydt:nrhs","Five inputs required.");
  }
  if(nlhs != n_outputs) {
    mexErrMsgIdAndTxt("bypass273_cvode_dydt:nlhs","One output required.");
  }


  if ( !mxIsDouble(prhs[input_y]) ) {
    mexErrMsgIdAndTxt("system_bypass273_cvode:not_double",
                      "Error: Input vector y is not a double.");
  }

  if ( mxGetNumberOfElements(prhs[input_y]) != system_size ) {
    mexErrMsgIdAndTxt("system_bypass273_cvode:wrong_size",
                     "Input vector y must have 273 elements.");
  }

  for ( int i = 0; i < n_parameters; i++ ) {
    if ( !mxIsDouble(prhs[i + 2]) ) {
      mexErrMsgIdAndTxt("bypass273_cvode_dydt:not_double",
                        "Error: one of the parameters is not a double");
    } else if ( mxIsComplex(prhs[i + 2]) ) {
      mexErrMsgIdAndTxt("bypass273_cvode_dydt:not_complex",
                        "Error: one of the parameters not a real number");
    } else if ( mxGetNumberOfElements(prhs[i + 2]) != 1 ) {
      mexErrMsgIdAndTxt("bypass273_cvode_dydt:notScalar",
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
  plhs[0] = mxCreateDoubleMatrix(system_size, n_parameters, mxREAL);

  /* get a pointer to the real data in the output matrix */
  #if MX_HAS_INTERLEAVED_COMPLEX
  outMatrix = mxGetDoubles(plhs[0]);
  #else
  outMatrix = mxGetPr(plhs[0]);
  #endif

  /* call the computational routine */
  bypass273_cvode_jacobian_params(y,parameters,outMatrix);
    
}

#define JAC(i,j) jac[(i) + (j) * system_size]

void bypass273_cvode_jacobian_params(double* y, double* parameters, double* jac) {
  JAC(268,0) =  y[269]*(3.9E+1/2.0E+2)+7.15709E-1;
  JAC(269,0) =  y[269]*(-3.9E+1/2.0E+2);
}