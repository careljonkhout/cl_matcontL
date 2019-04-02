/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * _coder_gauss_elim_Jacobian_fun_api.h
 *
 * Code generation for function '_coder_gauss_elim_Jacobian_fun_api'
 *
 */

#ifndef _CODER_GAUSS_ELIM_JACOBIAN_FUN_API_H
#define _CODER_GAUSS_ELIM_JACOBIAN_FUN_API_H

/* Include files */
#include "tmwtypes.h"
#include "mex.h"
#include "emlrt.h"
#include <stddef.h>
#include <stdlib.h>
#include "_coder_gauss_elim_Jacobian_fun_api.h"

/* Type Definitions */
#ifndef struct_emxArray_real_T
#define struct_emxArray_real_T

struct emxArray_real_T
{
  real_T *data;
  int32_T *size;
  int32_T allocatedSize;
  int32_T numDimensions;
  boolean_T canFreeData;
};

#endif                                 /*struct_emxArray_real_T*/

#ifndef typedef_emxArray_real_T
#define typedef_emxArray_real_T

typedef struct emxArray_real_T emxArray_real_T;

#endif                                 /*typedef_emxArray_real_T*/

/* Variable Declarations */
extern emlrtCTX emlrtRootTLSGlobal;
extern emlrtContext emlrtContextGlobal;

/* Function Declarations */
extern void gauss_elim_Jacobian_fun(emxArray_real_T *blocks, real_T ntst, real_T
  ncol, real_T nphase);
extern void gauss_elim_Jacobian_fun_api(const mxArray * const prhs[4], int32_T
  nlhs);
extern void gauss_elim_Jacobian_fun_atexit(void);
extern void gauss_elim_Jacobian_fun_initialize(void);
extern void gauss_elim_Jacobian_fun_terminate(void);
extern void gauss_elim_Jacobian_fun_xil_terminate(void);

#endif

/* End of code generation (_coder_gauss_elim_Jacobian_fun_api.h) */
