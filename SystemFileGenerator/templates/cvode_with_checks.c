#include <cvodes/cvodes.h>
#include <sunmatrix/sunmatrix_band.h>  /* access to band SUNMatrix       */
#include <sunlinsol/sunlinsol_band.h>  /* access to band SUNLinearSolver */
#include <sunmatrix/sunmatrix_dense.h>  /* access to dense SUNMatrix       */
#include <sunlinsol/sunlinsol_dense.h>  /* access to dense SUNLinearSolver */

static void check_null(void* ptr, const char* funcname) {
  if (ptr == NULL) {
    mexErrMsgIdAndTxt("CVODE:null_pointer",
             "SUNDIALS_ERROR: %s() failed - returned NULL pointer\n", funcname); 
  }
}

static void check(int flag, const char* funcname) {
  if (flag < 0) {
    mexErrMsgIdAndTxt("CVODE:error", 
                      "SUNDIALS_ERROR: %s() failed with error_code %s",
                      funcname, CVodeGetReturnFlagName(flag));
  }
}

void* CVodeCreate_nullchecked(int linear_multistep_method) {
  void* cvode = CVodeCreate(linear_multistep_method);
  check_null(cvode, "CVodeCreate_nullchecked");
  return cvode;
}

void CVodeSetErrHandlerFn_checked(void*          cvode, 
                                  CVErrHandlerFn f,
                                  void*          eh_data) {
  int flag = CVodeSetErrHandlerFn(cvode, f, eh_data);
  check(flag, "CVodeSetErrHandlerFn_checked");
}

void CVodeInit_checked(void* cvode, CVRhsFn rhs, realtype t0, N_Vector y0) {
  int flag = CVodeInit(cvode, rhs, t0, y0);
  check(flag, "CVodeInit_checked");
}

void CVodeSetUserData_checked(void* cvode, void* user_data) {
  int flag = CVodeSetUserData(cvode, user_data);
  check(flag, "CVodeSetUserData_checked");
}

void CVodeSStolerances_checked(void* cvode, realtype reltol, realtype abstol) {
  int flag = CVodeSStolerances(cvode, reltol, abstol);
  check(flag, "CVodeSStolerances_checked");
}

SUNMatrix SUNDenseMatrix_nullchecked(sunindextype m, sunindextype n) {
  SUNMatrix mat = SUNDenseMatrix(m, n);
  check_null(mat, "SUNDenseMatrix_nullchecked");
  return mat;
}

SUNLinearSolver SUNLinSol_Dense_nullchecked(N_Vector y, SUNMatrix m) {
  SUNLinearSolver s = SUNLinSol_Dense(y, m);
  check_null(s, "SUNLinSol_Dense_nullchecked");
  return s;
}

SUNMatrix SUNBandMatrix_nullchecked(sunindextype n, 
                                    sunindextype mu,
                                    sunindextype ml) {
  SUNMatrix m = SUNBandMatrix(n, mu, ml);
  check_null(m, "SUNBandMatrix_nullchecked");
  return m;
}

SUNLinearSolver SUNLinSol_Band_nullchecked(N_Vector y, SUNMatrix m) {
  SUNLinearSolver s = SUNLinSol_Band(y, m);
  check_null(s, "SUNLinSol_Band_nullchecked");
  return s;
}

void CVodeSetLinearSolver_checked(void* cvode, SUNLinearSolver s, SUNMatrix m) {
  int flag = CVodeSetLinearSolver(cvode, s, m);
  check(flag, "CVodeSetLinearSolver_checked");
}

void CVodeSetJacFn_checked(void* cvode, CVLsJacFn j) {
  int flag = CVodeSetJacFn(cvode, j);
  check(flag, "CVodeSetJacFn_checked");
}

void CVodeSetMaxNumSteps_checked(void* cvode, long int mxsteps) {
  int flag = CVodeSetMaxNumSteps(cvode, mxsteps);
  check(flag, "CVodeSetMaxNumSteps_checked");
}