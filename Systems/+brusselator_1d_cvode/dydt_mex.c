// this mex function computes the right hand side of system of ODEs called
//                       "brusselator_1d_cvode"

#include <math.h>
#include "mex.h"
#include "user_data.h"

#define N_INPUTS     N_PARAMETERS + 2

#define INPUT_T 0
#define INPUT_Y 1

void dydt_ode(double* y, double* dydt, double* parameters);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

  double parameters[N_PARAMETERS];  
  double *y;                      /* 1xN input matrix */
  size_t ncols;                   /* size of matrix */
  double *outMatrix;              /* output matrix */



  for ( int i = 0; i < N_PARAMETERS; i++ ) {
    if ( !mxIsDouble(prhs[i + 2]) ) {
      mexErrMsgIdAndTxt("dydt_mex:not_double",
                        "Error: one of the parameters is not a double");
    } else if ( mxIsComplex(prhs[i + 2]) ) {
      mexErrMsgIdAndTxt("dydt_mex:not_complex",
                        "Error: one of the parameters not a real number");
    } else if ( mxGetNumberOfElements(prhs[i + 2]) != 1 ) {
      mexErrMsgIdAndTxt("dydt_mex:notScalar",
                        "Error: one of the parameters is not scalar");
    } else {
      parameters[i] = mxGetScalar(prhs[i+2]);
    }
  }
  
  int neq = PDE_DIMENSION * ((int) parameters[0]);
    /* check for proper number of arguments */
  if(nrhs != N_INPUTS) {
    mexErrMsgIdAndTxt("dydt_mex:nrhs", "%d inputs required.", N_INPUTS);
  }

  if ( !mxIsDouble(prhs[INPUT_Y]) ) {
    mexErrMsgIdAndTxt("dydt_mex:not_double",
                      "Error: Input vector y is not a double.");
  }

  if ( mxGetNumberOfElements(prhs[INPUT_Y]) != neq ) {
    mexErrMsgIdAndTxt("dydt_mex:wrong_size",
              "Input vector y must have %d elements, but has %d elements",
              neq, mxGetNumberOfElements(prhs[INPUT_Y]));
  }
    
    
  /* get a pointer to the data in the input matrix  */
  #if MX_HAS_INTERLEAVED_COMPLEX
  y = mxGetDoubles(prhs[INPUT_Y]);
  #else
  y = mxGetPr(prhs[INPUT_Y]);
  #endif

  /* create the output matrix */
  plhs[0] = mxCreateDoubleMatrix(neq, 1, mxREAL);

  /* get a pointer to the data in the output matrix */
  #if MX_HAS_INTERLEAVED_COMPLEX
  outMatrix = mxGetDoubles(plhs[0]);
  #else
  outMatrix = mxGetPr(plhs[0]);
  #endif

  /* call the routine that computes the values of dydt */
  dydt_ode(y,outMatrix,parameters);
    
}
