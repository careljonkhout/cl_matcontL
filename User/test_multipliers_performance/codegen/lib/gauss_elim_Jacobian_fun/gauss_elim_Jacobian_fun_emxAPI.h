/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * gauss_elim_Jacobian_fun_emxAPI.h
 *
 * Code generation for function 'gauss_elim_Jacobian_fun_emxAPI'
 *
 */

#ifndef GAUSS_ELIM_JACOBIAN_FUN_EMXAPI_H
#define GAUSS_ELIM_JACOBIAN_FUN_EMXAPI_H

/* Include files */
#include <stddef.h>
#include <stdlib.h>
#include "rtwtypes.h"
#include "gauss_elim_Jacobian_fun_types.h"

/* Function Declarations */
extern emxArray_real_T *emxCreateND_real_T(int numDimensions, int *size);
extern emxArray_real_T *emxCreateWrapperND_real_T(double *data, int
  numDimensions, int *size);
extern emxArray_real_T *emxCreateWrapper_real_T(double *data, int rows, int cols);
extern emxArray_real_T *emxCreate_real_T(int rows, int cols);
extern void emxDestroyArray_real_T(emxArray_real_T *emxArray);
extern void emxInitArray_real_T(emxArray_real_T **pEmxArray, int numDimensions);

#endif

/* End of code generation (gauss_elim_Jacobian_fun_emxAPI.h) */
