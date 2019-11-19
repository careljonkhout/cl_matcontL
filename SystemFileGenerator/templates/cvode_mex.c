#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <mex.h>

#include <cvodes/cvodes.h>
#include <nvector/nvector_serial.h>
//#include <nvector/nvector_openmp.h>
#include <sunmatrix/sunmatrix_band.h>   /* access to band SUNMatrix       */
#include <sunlinsol/sunlinsol_band.h>   /* access to band SUNLinearSolver */
#include <sunmatrix/sunmatrix_dense.h>  /* access to dense SUNMatrix       */
#include <sunlinsol/sunlinsol_dense.h>  /* access to dense SUNLinearSolver */

#include <sundials/sundials_types.h>
#include <sundials/sundials_math.h>

#include "user_data.h"
#include "cvode_with_checks.c"
#include "mex_utils.c"

#define OUTPUT_T           0
#define OUTPUT_Y           1
#define OUTPUT_SENSITIVITY 2
#define DEFAULT_ABS_TOL    RCONST(1.e-6) /* default scalar absolute tolerance */
#define DEFAULT_REL_TOL    RCONST(1.e-6) /* default relative tolerance */
#define NS                 1
#define INCREASING         1
#define MAX_NUM_STEPS      10*1000*1000

// implemented in dydt_cvode.c
int dydt_cvode(realtype t, N_Vector u, N_Vector udot, void *);

#if ANALYTIC_JACOBIAN
// implemented in jacobian_cvode.c
int jacobian_dydt(
               realtype t, N_Vector y_vec, N_Vector fy, SUNMatrix jac_structure,
              void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
#endif
// implemented in d_sensitivity_dt.c
int d_sensitivity_dt(int Ns, realtype t,
                N_Vector y_vector, 
                N_Vector ydot,
                int iS,
                N_Vector s_vector, 
                N_Vector ds_vector,
                void *user_data,
                N_Vector tmp1, N_Vector tmp2);


int return_to_plane(realtype t, N_Vector u, realtype* g  , void *);

void error_handler(int error_code, const char *module,
                          const char *function, char *msg, void *eh_data);

N_Vector get_N_Vector(const mxArray* array, char* input_name, int size);
N_Vector new_N_Vector(int length);

void check_argument_given(void* pointer,                char* arrayname);

void PrintFinalStats(void *cvode, booleantype sensi,
                            booleantype err_con, int sensi_meth);


realtype* get_reals(const mxArray* array, char* input_name, int size);
void copy_real_to_double(realtype* r, double* d, int length);
void mex_eval(const char * format, ... );


realtype     lower_bound_period;
N_Vector     point_on_limitcycle;
N_Vector     tangent_to_limitcycle;

void mexFunction(int n_output,       mxArray *mex_output[], 
                 int n_input,  const mxArray *mex_input []  ) {
  
  int          n_points                 = 0;
  realtype*    t_values                 = NULL;
  realtype*    initial_point            = NULL;

  int          sensitivity              = false;
  realtype*    sensitivity_vector       = NULL;

  booleantype  cycle_detection          = false;
               lower_bound_period       = 0;
  booleantype  lower_bound_period_given = false;
               point_on_limitcycle      = NULL;
               tangent_to_limitcycle    = NULL; 

  booleantype  verbose                  = false;
  
  realtype     abs_tol                  = DEFAULT_ABS_TOL;
  realtype     rel_tol                  = DEFAULT_REL_TOL;
  /* abs_tol and rel_tol will be overwritten if the caller specifies them */
  
  realtype     initial_step_size        = 0;
  booleantype  initial_step_size_given  = false;
  
  
  UserData user_data = malloc(sizeof(UserData));
  check_null(user_data, "malloc user_data");
  user_data->parameters = NULL;
  

  
  if (n_input % 2 != 0) {
    mexErrMsgIdAndTxt("cvode:nrhs_uneven",
                     "The number of arguments is uneven. "
                     "The number for arguments is %d. "
                     "Arguments must be supplied as name-value pairs", n_input);
  }

  for ( int i = 0; i < n_input; i += 2 ) {
    
    if ( !mxIsChar(mex_input[i]) ) {
      mexErrMsgIdAndTxt("cvode:not_name_value", 
                              "arguments must be supplied in name-value pairs");
    }
    char*           arg_name = (char*) mxArrayToString(mex_input[i]);
    const mxArray*  value    = mex_input[i + 1];

    if        ( ! strcmp(arg_name, "t_values"             ) ) {
      
      
      n_points = mxGetNumberOfElements(value);
      t_values = get_reals(value, "t_values", n_points);
      
      if (n_points < 2) {
        mexErrMsgIdAndTxt("cvode:wrong_size",
                        "Input vector t_values must have at least 2 elements.");
      }
      for (int i = 0; i < n_points - 1; i ++) {
        if ( t_values[i] >= t_values[i+1] ) {
          mexErrMsgIdAndTxt("cvode:not_increasing",
                 "The sequence of time points must be an increasing sequence.");
        }
      }
      
    } else if ( ! strcmp(arg_name, "initial_point"        ) ) {
      
      initial_point         = get_reals(value, "initial_point", NEQ);
      
    } else if ( ! strcmp(arg_name, "ode_parameters"       ) ) {
            
      user_data->parameters = get_reals(value, "ode_parameters", N_PARAMETERS);
      
    } else if ( ! strcmp(arg_name, "sensitivity_vector"   ) ) {
      
      sensitivity_vector    = get_reals(value, "sensitivity_vector", NEQ);
      
      sensitivity = true;
      
    } else if ( ! strcmp(arg_name, "cycle_detection"      ) ) {
      
      check_bool(value, "cycle_detection");
      cycle_detection = mxIsLogicalScalarTrue(value);
      
    } else if ( ! strcmp(arg_name, "lower_bound_period"  ) ) {
      
      lower_bound_period = get_positive_scalar(value, "lower_bound_period");
      lower_bound_period_given = true;
      
    } else if ( ! strcmp(arg_name, "point_on_limitcycle"  ) ) {
      
      point_on_limitcycle = get_N_Vector(value, "point_on_limitcycle", NEQ);
      tangent_to_limitcycle = new_N_Vector(NEQ);
      
    } else if ( ! strcmp(arg_name, "abs_tol"  ) ) {
      
      abs_tol = get_positive_scalar(value, "abs_tol");
      
    } else if ( ! strcmp(arg_name, "rel_tol"  ) ) {
      
      rel_tol = get_positive_scalar(value, "rel_tol");
      
    } else if ( ! strcmp(arg_name, "initial_step_size"  ) ) {
      
      initial_step_size = get_positive_scalar(value, "initial step size");
      initial_step_size_given = true;
      
    } else if ( ! strcmp(arg_name, "verbose"      ) ) {
      
      check_bool(value, "verbose");
      verbose = mxIsLogicalScalarTrue(value);
      
    } else {
  
      mexErrMsgIdAndTxt("cvode:invalid_argument_name", 
                        "Invalid argument name for cvode: %s.", arg_name);
    }
  }
  

  if(n_output <  2) {
    mexErrMsgIdAndTxt("cvode:n_output", "2 outputs are needed: t and y "
                     "(or 3 outputs for sensitivity analysis)");
  }
  
  check_argument_given((void *) t_values               , "t_values"     );
  check_argument_given((void *) initial_point          , "initial_point");
  check_argument_given((void *) (user_data->parameters), "parameters"   );
  

  if (cycle_detection) {
    // we wait to compute tangent_to_limitcycle, since we need parameters to 
    // be loaded first
    dydt_cvode(0, point_on_limitcycle, tangent_to_limitcycle, user_data);
    
    if (sensitivity) {
       mexErrMsgIdAndTxt("cvode:sens_and_cd",
                         "simultaneous cycle detection and "
                         "sensitivity analysis is not supported");   
    }
    if (! lower_bound_period_given) {
      mexErrMsgIdAndTxt("cvode:no_lower_bound_period",
                        "You enabled cycle detection, "
                        "but you did not specify lower_bound_period");
    }
    if (point_on_limitcycle == NULL) {
      mexErrMsgIdAndTxt("cvode:no_point_on_limitcycle",
                        "You enabled cycle detection, "
                        "but you did not specify point_on_limit_cycle");    
    }
  }
  
  if(sensitivity && n_output != 3) {
    mexErrMsgIdAndTxt("cvode:n_output",
                      "You enabled sensitivity analysis, " 
                      "but you specified less than three outputs. "
                      "Three outputs are needed for sensitivity analysis: "
                      "t (time), y (state variables), and s (sensitivity");
  }
  

  void*             cvode       = NULL;
  SUNMatrix         A           = NULL;
  SUNLinearSolver   LS          = NULL;
  int               sensi_meth  = CV_SIMULTANEOUS;
  booleantype       err_con     = true;
    
  N_Vector solver_y = new_N_Vector(NEQ);
  N_Vector solver_s;
  if (sensitivity) {
    solver_s = new_N_Vector(NEQ);
  }

  
  realtype* y_data             = N_VGetArrayPointer(solver_y);
  
  memcpy(y_data, initial_point, NEQ * sizeof(realtype));
  
  realtype* s_data;
  if ( sensitivity ) {
    s_data = N_VGetArrayPointer(solver_s); 
    memcpy(s_data, sensitivity_vector, NEQ * sizeof(realtype));
  }
  
  cvode = CVodeCreate_nullchecked(CV_BDF);

  CVodeInit_checked           (cvode, dydt_cvode,    t_values[0], solver_y);  
  CVodeSetErrHandlerFn_checked(cvode, error_handler, NULL                 );
  CVodeSetUserData_checked    (cvode, user_data                           );
  CVodeSStolerances_checked   (cvode, rel_tol,       abs_tol              );
  CVodeSetMaxNumSteps_checked (cvode, MAX_NUM_STEPS                       );
  CVodeSetMaxHnilWarns        (cvode, 1                                   );
  #ifdef MAX_ORDER_BDF
  CVodeSetMaxOrd              (cvode, MAX_ORDER_BDF);
  #endif
  #ifdef SUNDIALS_EXTENDED_PRECISION
  CVodeSetMinStep             (cvode, 1e-20q); // q means quadruple precision
  #endif
  
  if (initial_step_size_given) {
    CVodeSetInitStep(cvode, initial_step_size);
  }
     
  #if JACOBIAN_STORAGE == DENSE
  A  = SUNDenseMatrix_nullchecked(NEQ, NEQ);
  LS = SUNLinSol_Dense_nullchecked(solver_y, A);
  #endif
  
  #if JACOBIAN_STORAGE == BANDED
  A  = SUNBandMatrix_nullchecked(NEQ, UPPER_BANDWIDTH, LOWER_BANDWIDTH);
  LS = SUNLinSol_Band_nullchecked(solver_y, A);
  #endif
   
  CVodeSetLinearSolver_checked(cvode, LS, A);
 
  /* Set the user-supplied Jacobian routine Jac */
  #if ANALYTIC_JACOBIAN
  CVodeSetJacFn(cvode, jacobian_dydt);
  #endif
  
  int rootfinding_direction = INCREASING;
  if (cycle_detection) {
    CVodeRootInit        (cvode, 1,                     return_to_plane);
    CVodeSetRootDirection(cvode, &rootfinding_direction                );
  }
  
  if (sensitivity) {
    #if ANALYTIC_JACOBIAN
    CVodeSensInit1       (cvode, NS, sensi_meth, d_sensitivity_dt, &solver_s);
    #else
    CVodeSensInit1       (cvode, NS, sensi_meth, NULL            , &solver_s);
    #endif
    CVodeSensEEtolerances(cvode);
    //CVodeSensSStolerances(cvode, rel_tol, &abs_tol);
    CVodeSetSensErrCon   (cvode, true);    
    realtype dummy_value = 0;
    CVodeSetSensParams(cvode, &dummy_value, NULL, NULL);
    if (verbose) {
      mexPrintf("Sensitivity: YES ");
      switch (sensi_meth) {
        case CV_SIMULTANEOUS: mexPrintf("( SIMULTANEOUS +"); break;
        case CV_STAGGERED   : mexPrintf("( STAGGERED +"   ); break;
        case CV_STAGGERED1  : mexPrintf("( STAGGERED1 +"  ); break;
      }
      mexPrintf(err_con ? " FULL ERROR CONTROL )" : " PARTIAL ERROR CONTROL )");
      mexPrintf("\n");
    }
  } else {
    if (verbose) mexPrintf("Sensitivity: NO \n");
  }
  
  
  
  mxArray* t_output_array = mxCreateDoubleMatrix(1,   n_points, mxREAL);
  mxArray* y_output_array = mxCreateDoubleMatrix(NEQ, n_points, mxREAL);
  if (sensitivity) {
    mex_output[OUTPUT_SENSITIVITY] = mxCreateDoubleMatrix(NEQ, 1, mxREAL);
  }
  
  mxDouble* t_out = my_mex_get_doubles(t_output_array);
  mxDouble* y_out = my_mex_get_doubles(y_output_array);
  mxDouble* s_out;
  if (sensitivity) {
    s_out = my_mex_get_doubles(mex_output[OUTPUT_SENSITIVITY]);
  }
  
  copy_real_to_double(t_values, t_out,  n_points);
  copy_real_to_double(y_data,   y_out, NEQ      );
  
  double* y_out_ptr = y_out + NEQ;



  realtype t          = 0;
  int      rootsfound = false;
  int      i;
  
  for(i = 1; i < n_points; i++) {
    // mex_eval("print_temp('t=%.15e')",(double)t_values[i]);
    int flag = CVode(cvode, t_values[i], solver_y, &t, CV_NORMAL);
    check(flag, "CVode");
    if (sensitivity) {
      flag = CVodeGetSens(cvode, &t, &solver_s);
      check(flag, "CVodeGetSens");
      copy_real_to_double(s_data, s_out, NEQ);
      // todo: fix: s_out seems to be overwritten n_points times here
    }
    if (cycle_detection) {
      CVodeGetRootInfo(cvode, &rootsfound);
    }
    copy_real_to_double(y_data, y_out_ptr, NEQ);
    y_out_ptr += NEQ;

    if (rootsfound) {
      //mexPrintf("cvode: period: %.6f\n", t);
      t_out[i] = t;
      break;
    }
  }
  
  if (i+1 < n_points) {
    // this happens when a cycle has been detected
    mxSetN(t_output_array, i+1);
    mxSetN(y_output_array, i+1);
  }
  mxArray* y_transposed = mxCreateDoubleMatrix(i+1, NEQ, mxREAL);
  mexCallMATLAB(1, &y_transposed, 1, &y_output_array, "transpose");
  mex_output[OUTPUT_Y] = y_transposed;
  
  mxArray* t_transposed = mxCreateDoubleMatrix(i+1, 1, mxREAL);
  mexCallMATLAB(1, &t_transposed, 1, &t_output_array, "transpose");
  mex_output[OUTPUT_T] = t_transposed;
  
   
  if (verbose) {
    PrintFinalStats(cvode, sensitivity, err_con, sensi_meth);
  }
  
  /* Free working memory. Note that matlab would take care of it automatically 
   * if we do not do it here, but freeing memory as soon as we can is better. */
  
  mxDestroyArray(y_output_array);
  mxDestroyArray(t_output_array);

  N_VDestroy(solver_y);
  if (sensitivity) {
    N_VDestroy(solver_s);
  }
  if (cycle_detection) {
    N_VDestroy(point_on_limitcycle);
    N_VDestroy(tangent_to_limitcycle);
  }
  
  SUNLinSolFree(LS); /* Free the linear solver memory */
  SUNMatDestroy(A);  /* Free the matrix memory */
  free(user_data);   /* Free the user data */

  CVodeFree(&cvode);
}

int return_to_plane(realtype t, N_Vector y, realtype *gout, void *user_data) {
  if (t > lower_bound_period) {
    *gout = 0;
    realtype* y_data   = NV_DATA_S(y);
    realtype* pol_data = NV_DATA_S(point_on_limitcycle);
    realtype* ttl_data = NV_DATA_S(tangent_to_limitcycle);
    for (int i = 0; i < NV_LENGTH_S(y); i++) {
      *gout += (y_data[i] - pol_data[i]) * ttl_data[i];
    }
  } else {
    *gout = 1;
  }
  return 0;
}

void error_handler(int error_code, const char* module,
         const char* function, char* message, void *eh_data) {
  if (error_code < 0) {
    mexErrMsgIdAndTxt("cvode:integrator_error",
            "Error in %s() in sundails module %s. "
            "error code: %d: %s.\n  %s\n",
            function, module,
            error_code, CVodeGetReturnFlagName(error_code), message);
  } else {
    mex_eval("print_diag(1, '%s\\n')", message); ;
  }
}

void check_argument_given(void * pointer, char* argument_name) {
  if (argument_name == NULL) {
    mexErrMsgIdAndTxt("cvode:missing_argument",  
                      "Argument '%s' was not passed to cvode.", argument_name);  
  }
}

N_Vector new_N_Vector(int length) {
  N_Vector v;
  #if OPEN_MP
  v = N_VNew_OpenMP(NEQ, N_THREADS); 
  #else
  v = N_VNew_Serial(NEQ); 
  #endif
  return v;
}

N_Vector get_N_Vector(const mxArray* array, char* input_name, int size) {
  realtype* data = get_reals(array, input_name, size);
  return N_VMake_Serial(size, data);
}



/*
 * Print some final statistics located in the CVODES memory
 */

void PrintFinalStats(void *cvode, booleantype sensi,
                            booleantype err_con, int sensi_meth) {
  long int nst;
  long int nfe, nsetups, nni, ncfn, netf;
  long int nfSe, nfeS, nsetupsS, nniS, ncfnS, netfS;
  int flag;

  flag = CVodeGetNumSteps(cvode, &nst);
  check(flag, "CVodeGetNumSteps");
  flag = CVodeGetNumRhsEvals(cvode, &nfe);
  check(flag, "CVodeGetNumRhsEvals");
  flag = CVodeGetNumLinSolvSetups(cvode, &nsetups);
  check(flag, "CVodeGetNumLinSolvSetups");
  flag = CVodeGetNumErrTestFails(cvode, &netf);
  check(flag, "CVodeGetNumErrTestFails");
  flag = CVodeGetNumNonlinSolvIters(cvode, &nni);
  check(flag, "CVodeGetNumNonlinSolvIters");
  flag = CVodeGetNumNonlinSolvConvFails(cvode, &ncfn);
  check(flag, "CVodeGetNumNonlinSolvConvFails");

  if (sensi) {
    flag = CVodeGetSensNumRhsEvals(cvode, &nfSe);
    check(flag, "CVodeGetSensNumRhsEvals");
    flag = CVodeGetNumRhsEvalsSens(cvode, &nfeS);
    check(flag, "CVodeGetNumRhsEvalsSens");
    flag = CVodeGetSensNumLinSolvSetups(cvode, &nsetupsS);
    check(flag, "CVodeGetSensNumLinSolvSetups");
    if (err_con) {
      flag = CVodeGetSensNumErrTestFails(cvode, &netfS);
      check(flag, "CVodeGetSensNumErrTestFails");
    } else {
      netfS = 0;
    }
    if ((sensi_meth == CV_STAGGERED) || (sensi_meth == CV_STAGGERED1)) {
      flag = CVodeGetSensNumNonlinSolvIters(cvode, &nniS);
      check(flag, "CVodeGetSensNumNonlinSolvIters");
      flag = CVodeGetSensNumNonlinSolvConvFails(cvode, &ncfnS);
      check(flag, "CVodeGetSensNumNonlinSolvConvFails");
    } else {
      nniS = 0;
      ncfnS = 0;
    }
  }

  mexPrintf("CVode Statistics: ");
  mexPrintf("number of steps = %d ", nst);
  mexPrintf("number of rhs evals = %d\n",   nfe);
  //mexPrintf("netf    = %5ld    nsetups  = %d\n", netf, nsetups);
  //mexPrintf("nni     = %5ld    ncfn     = %5\n", nni, ncfn);

  if(sensi) {
    mexPrintf("\n");
    mexPrintf("nfSe    = %5ld    nfeS     = %5ld\n", nfSe, nfeS);
    mexPrintf("netfs   = %5ld    nsetupsS = %5ld\n", netfS, nsetupsS);
    mexPrintf("nniS    = %5ld    ncfnS    = %5ld\n", nniS, ncfnS);
  }

}

realtype* get_reals(const mxArray* array, char* input_name, int size) {
  check_double(array,       input_name);
  check_real  (array,       input_name);
  check_size  (array, size, input_name);
  
  #ifdef SUNDIALS_EXTENDED_PRECISION
  realtype* real_data   = malloc(size * sizeof(realtype));
  double*   matlab_data = my_mex_get_doubles(array);
  for (int i = 0; i < size; i++) {
    real_data[i] = (realtype) matlab_data[i];
  }
  return real_data;
  #else
  double* matlab_data = my_mex_get_doubles(array);
  return matlab_data;
  #endif
}

void copy_real_to_double(realtype* r, double* d, int length) {
  #ifdef SUNDIALS_EXTENDED_PRECISION
  for (int i = 0; i < length; i++) {
    d[i] = (double) r[i];
  }
  #else
  memcpy(d, r, length * sizeof(double));
  #endif
}

void mex_eval(const char * format, ... ) {
  va_list args1;
  va_start(args1, format);
  va_list args2;
  va_copy(args2, args1);
  char buffer[1 + vsnprintf(NULL, 0, format, args1)];
  va_end(args1);
  vsnprintf(buffer, sizeof buffer, format, args2);
  va_end(args2);
  mexEvalString(buffer);
}