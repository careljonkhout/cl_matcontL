/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * main.c
 *
 * Code generation for function 'main'
 *
 */

/*************************************************************************/
/* This automatically generated example C main file shows how to call    */
/* entry-point functions that MATLAB Coder generated. You must customize */
/* this file for your application. Do not modify this file directly.     */
/* Instead, make a copy of this file, modify it, and integrate it into   */
/* your development environment.                                         */
/*                                                                       */
/* This file initializes entry-point function arguments to a default     */
/* size and value before calling the entry-point functions. It does      */
/* not store or use any values returned from the entry-point functions.  */
/* If necessary, it does pre-allocate memory for returned values.        */
/* You can use this file as a starting point for a main function that    */
/* you can deploy in your application.                                   */
/*                                                                       */
/* After you copy the file, and before you deploy it, you must make the  */
/* following changes:                                                    */
/* * For variable-size function arguments, change the example sizes to   */
/* the sizes that your application requires.                             */
/* * Change the example values of function arguments to the values that  */
/* your application requires.                                            */
/* * If the entry-point functions return values, store these values or   */
/* otherwise use them as required by your application.                   */
/*                                                                       */
/*************************************************************************/
/* Include files */
#include "gauss_elim_Jacobian_fun.h"
#include "main.h"
#include "gauss_elim_Jacobian_fun_terminate.h"
#include "gauss_elim_Jacobian_fun_emxAPI.h"
#include "gauss_elim_Jacobian_fun_initialize.h"

/* Function Declarations */
static double argInit_real_T(void);
static emxArray_real_T *c_argInit_UnboundedxUnboundedxU(void);
static void main_gauss_elim_Jacobian_fun(void);

/* Function Definitions */
static double argInit_real_T(void)
{
  return 0.0;
}

static emxArray_real_T *c_argInit_UnboundedxUnboundedxU(void)
{
  emxArray_real_T *result;
  static int iv0[3] = { 2, 2, 2 };

  int idx0;
  int idx1;
  int idx2;

  /* Set the size of the array.
     Change this size to the value that the application requires. */
  result = emxCreateND_real_T(3, iv0);

  /* Loop over the array to initialize each element. */
  for (idx0 = 0; idx0 < result->size[0U]; idx0++) {
    for (idx1 = 0; idx1 < result->size[1U]; idx1++) {
      for (idx2 = 0; idx2 < result->size[2U]; idx2++) {
        /* Set the value of the array element.
           Change this value to the value that the application requires. */
        result->data[(idx0 + result->size[0] * idx1) + result->size[0] *
          result->size[1] * idx2] = argInit_real_T();
      }
    }
  }

  return result;
}

static void main_gauss_elim_Jacobian_fun(void)
{
  emxArray_real_T *blocks;

  /* Initialize function 'gauss_elim_Jacobian_fun' input arguments. */
  /* Initialize function input argument 'blocks'. */
  blocks = c_argInit_UnboundedxUnboundedxU();

  /* Call the entry-point 'gauss_elim_Jacobian_fun'. */
  gauss_elim_Jacobian_fun(blocks, argInit_real_T(), argInit_real_T(),
    argInit_real_T());
  emxDestroyArray_real_T(blocks);
}

int main(int argc, const char * const argv[])
{
  (void)argc;
  (void)argv;

  /* Initialize the application.
     You do not need to do this more than one time. */
  gauss_elim_Jacobian_fun_initialize();

  /* Invoke the entry-point functions.
     You can call entry-point functions multiple times. */
  main_gauss_elim_Jacobian_fun();

  /* Terminate the application.
     You do not need to do this more than one time. */
  gauss_elim_Jacobian_fun_terminate();
  return 0;
}

/* End of code generation (main.c) */
