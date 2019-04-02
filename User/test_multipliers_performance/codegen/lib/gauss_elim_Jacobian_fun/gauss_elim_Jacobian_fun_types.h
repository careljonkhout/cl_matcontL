/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * gauss_elim_Jacobian_fun_types.h
 *
 * Code generation for function 'gauss_elim_Jacobian_fun'
 *
 */

#ifndef GAUSS_ELIM_JACOBIAN_FUN_TYPES_H
#define GAUSS_ELIM_JACOBIAN_FUN_TYPES_H

/* Include files */
#include "rtwtypes.h"

/* Type Definitions */
#ifndef struct_emxArray_real_T
#define struct_emxArray_real_T

struct emxArray_real_T
{
  double *data;
  int *size;
  int allocatedSize;
  int numDimensions;
  boolean_T canFreeData;
};

#endif                                 /*struct_emxArray_real_T*/

#ifndef typedef_emxArray_real_T
#define typedef_emxArray_real_T

typedef struct emxArray_real_T emxArray_real_T;

#endif                                 /*typedef_emxArray_real_T*/
#endif

/* End of code generation (gauss_elim_Jacobian_fun_types.h) */
