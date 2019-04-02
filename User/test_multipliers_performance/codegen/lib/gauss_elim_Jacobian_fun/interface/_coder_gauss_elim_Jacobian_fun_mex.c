/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * _coder_gauss_elim_Jacobian_fun_mex.c
 *
 * Code generation for function '_coder_gauss_elim_Jacobian_fun_mex'
 *
 */

/* Include files */
#include "_coder_gauss_elim_Jacobian_fun_api.h"
#include "_coder_gauss_elim_Jacobian_fun_mex.h"

/* Function Declarations */
static void c_gauss_elim_Jacobian_fun_mexFu(int32_T nlhs, int32_T nrhs, const
  mxArray *prhs[4]);

/* Function Definitions */
static void c_gauss_elim_Jacobian_fun_mexFu(int32_T nlhs, int32_T nrhs, const
  mxArray *prhs[4])
{
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  st.tls = emlrtRootTLSGlobal;

  /* Check for proper number of arguments. */
  if (nrhs != 4) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:WrongNumberOfInputs", 5, 12, 4, 4,
                        23, "gauss_elim_Jacobian_fun");
  }

  if (nlhs > 0) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:TooManyOutputArguments", 3, 4, 23,
                        "gauss_elim_Jacobian_fun");
  }

  /* Call the function. */
  gauss_elim_Jacobian_fun_api(prhs, nlhs);
}

void mexFunction(int32_T nlhs, mxArray *plhs[], int32_T nrhs, const mxArray
                 *prhs[])
{
  (void)plhs;
  mexAtExit(gauss_elim_Jacobian_fun_atexit);

  /* Module initialization. */
  gauss_elim_Jacobian_fun_initialize();

  /* Dispatch the entry-point. */
  c_gauss_elim_Jacobian_fun_mexFu(nlhs, nrhs, prhs);

  /* Module termination. */
  gauss_elim_Jacobian_fun_terminate();
}

emlrtCTX mexFunctionCreateRootTLS(void)
{
  emlrtCreateRootTLS(&emlrtRootTLSGlobal, &emlrtContextGlobal, NULL, 1);
  return emlrtRootTLSGlobal;
}

/* End of code generation (_coder_gauss_elim_Jacobian_fun_mex.c) */
