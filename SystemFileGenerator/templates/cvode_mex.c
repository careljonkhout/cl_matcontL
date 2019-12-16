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
#include "mex_utils_cvode.c"

#define OUTPUT_T            0
#define OUTPUT_Y            1
#define OUTPUT_SENSITIVITY  2
#define INCREASING          1

#define DEFAULT_ABS_TOL     RCONST(1.e-6) /* default scalar absolute tolerance */
#define DEFAULT_REL_TOL     RCONST(1.e-6) /* default relative tolerance */
#define MAX_NUM_STEPS       10*1000*1000
#define SENSITIVITY_METHOD  CV_SIMULTANEOUS

// implemented in dydt_cvode.c
int dydt_cvode(realtype t, N_Vector u, N_Vector udot, void *);

#if ANALYTIC_JACOBIAN

// implemented in jacobian_cvode.c
int jacobian_dydt(
               realtype t, N_Vector y_vec, N_Vector fy, SUNMatrix jac_structure,
              void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

// implemented in d_sensitivity_dt.c
int d_sensitivity_dt(int Ns, realtype t,
                N_Vector y_vector, 
                N_Vector ydot,
                int iS,
                N_Vector s_vector, 
                N_Vector ds_vector,
                void *user_data,
                N_Vector tmp1, N_Vector tmp2);

// implemented in d_sensitivity_dt_pars.c
int d_sensitivity_dt_pars(int Ns, realtype t,
                N_Vector y_vector, 
                N_Vector ydot,
                int iS,
                N_Vector s_vector, 
                N_Vector ds_vector,
                void *user_data,
                N_Vector tmp1, N_Vector tmp2);



#endif

int return_to_plane(realtype t, N_Vector u, realtype* g  , void *);

void error_handler(int error_code, const char *module,
                          const char *function, char *msg, void *eh_data);

N_Vector get_N_Vector(Argument arg, int size);
N_Vector new_N_Vector(int length);

void check_argument_given(void* pointer,                char* arrayname);

void PrintFinalStats(void *cvode, booleantype sensitivity);


realtype* get_reals(Argument arg, int size);
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

  booleantype  sensitivity              = false;
  int          n_sensitivity            = 0;
  realtype*    sensitivity_vectors      = NULL;
  
  booleantype  parameter_sensitivity    = false;
  int          parameter_index          = 0;

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
  Argument arg;
  

  
  if (n_input % 2 != 0) {
    mexErrMsgIdAndTxt("cvode:nrhs_uneven",
                     "The number of arguments is uneven. "
                     "The number for arguments is %d. "
                     "Arguments must be supplied as name-value pairs", n_input);
  }

  for ( int i = 0; i < n_input; i += 2 ) {
    
    if ( !mxIsChar(mex_input[i]) ) {
      mexErrMsgIdAndTxt("cvode:not_name_value", 
                              "Arguments must be supplied in name-value pairs,"
                              " but argument %d is not a name", i+1);
    }
    arg.name  = (char*) mxArrayToString(mex_input[i]);
    arg.value = mex_input[i + 1];

    if        ( ! strcmp(arg.name, "t_values"             ) ) {
      
      n_points = mxGetNumberOfElements(arg.value);
      t_values = get_reals(arg, n_points);
      
      if (n_points < 2) {
        mexErrMsgIdAndTxt("cvode:wrong_size",
                        "Input vector t_values must have at least 2 elements.");
      }
      for (int i = 0; i < n_points - 1; i ++) {
        if ( t_values[i] >= t_values[i + 1] ) {
          mexErrMsgIdAndTxt("cvode:not_increasing",
                 "The sequence of time points must be an increasing sequence.");
        }
      }
      
    } else if ( ! strcmp(arg.name, "initial_point"        ) ) {
      
      initial_point         = get_reals(arg, NEQ);
      
    } else if ( ! strcmp(arg.name, "ode_parameters"       ) ) {
            
      user_data->parameters = get_reals(arg, N_PARAMETERS);
      
    } else if ( ! strcmp(arg.name, "sensitivity_vector"   ) ) {
      
      int n_elements      = mxGetNumberOfElements(arg.value);
      sensitivity_vectors = get_reals(arg, n_elements);
      n_sensitivity       = mxGetN(arg.value);
      sensitivity         = true;
      
    } else if ( ! strcmp(arg.name, "cycle_detection"      ) ) {
      
      cycle_detection = get_bool(arg);
      
    } else if ( ! strcmp(arg.name, "lower_bound_period"   ) ) {
      
      lower_bound_period = get_positive_scalar(arg);
      lower_bound_period_given = true;
      
    } else if ( ! strcmp(arg.name, "point_on_limitcycle"  ) ) {
      
      point_on_limitcycle = get_N_Vector(arg, NEQ);
      tangent_to_limitcycle = new_N_Vector(NEQ);
      
    } else if ( ! strcmp(arg.name, "abs_tol"              ) ) {
      
      abs_tol = get_positive_scalar(arg);
      
    } else if ( ! strcmp(arg.name, "rel_tol"              ) ) {
      
      rel_tol = get_positive_scalar(arg);
      
    } else if ( ! strcmp(arg.name, "initial_step_size"    ) ) {
      
      initial_step_size = get_positive_scalar(arg);
      initial_step_size_given = true;
      
    } else if ( ! strcmp(arg.name, "verbose"              ) ) {
      
      verbose = get_bool(arg);
    
    } else if ( ! strcmp(arg.name, "parameter_sensitivity") ) {
      
      parameter_index = get_nonnegative_int(arg);
      n_sensitivity = 1;
      parameter_sensitivity = true;
      
    } else {
  
      mexErrMsgIdAndTxt("cvode:invalid_argument_name", 
                        "Invalid argument name for cvode: %s.", arg.name);
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
  
  if( ( sensitivity || parameter_sensitivity) && n_output != 3) {
    mexErrMsgIdAndTxt("cvode:n_output",
                      "You enabled sensitivity analysis, " 
                      "but you specified less than three outputs. "
                      "Three outputs are needed for sensitivity analysis: "
                      "t (time), y (state variables), and s (sensitivity)");
  }
  

  void*             cvode       = NULL;
  SUNMatrix         A           = NULL;
  SUNLinearSolver   LS          = NULL;
    
  N_Vector solver_y = new_N_Vector(NEQ);
  N_Vector* solver_s;
  realtype* s_data;
  if (sensitivity || parameter_sensitivity) {
    solver_s = malloc(n_sensitivity * sizeof(N_Vector *));
    if (sensitivity) {
      s_data = sensitivity_vectors;
      // create N_Vectors whose data pointers point to the vectors in 
      // sensitivity_vectors
      for (int i = 0; i < n_sensitivity; i++) {
        solver_s[i] = N_VMake_Serial(NEQ, sensitivity_vectors + i * NEQ);
      }
    } else if (parameter_sensitivity) {
      // allocate new N_Vector
      solver_s[0] = N_VNew_Serial(NEQ);
      s_data = N_VGetArrayPointer(solver_s[0]);
      for (int i = 0; i < NEQ; i++) {
        s_data[i] = 0.0;
      }
    }
  }

  N_VSetArrayPointer(initial_point, solver_y);
  realtype* y_data = initial_point;
  

  
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
  
  if (sensitivity || parameter_sensitivity) {
    #if ANALYTIC_JACOBIAN
    if (sensitivity) {
      CVodeSensInit1(cvode, n_sensitivity, SENSITIVITY_METHOD, 
                                              d_sensitivity_dt, solver_s);
    } else if (parameter_sensitivity) {
      CVodeSensInit1(cvode,             1, SENSITIVITY_METHOD, 
                                              d_sensitivity_dt_pars, solver_s);
    }
    #else
    CVodeSensInit1(cvode, NS, SENSITIVITY_METHOD, NULL, solver_s);
    #endif
    
    if (parameter_sensitivity) {
      CVodeSetSensParams(cvode, NULL, NULL, &parameter_index);
    }

    CVodeSensEEtolerances(cvode);
    CVodeSetSensErrCon(cvode, true);    
  } 
  
  
  
  mxArray* t_output_array = mxCreateDoubleMatrix(1,   n_points, mxREAL);
  mxArray* y_output_array = mxCreateDoubleMatrix(NEQ, n_points, mxREAL);
  if (sensitivity) {
    mex_output[OUTPUT_SENSITIVITY] = mxCreateDoubleMatrix(NEQ, n_sensitivity, 
                                                               mxREAL);
  }
    
  if (parameter_sensitivity) {
    mex_output[OUTPUT_SENSITIVITY] = mxCreateDoubleMatrix(NEQ, 1, mxREAL);
  }
  
  mxDouble* t_out = my_mex_get_doubles(t_output_array);
  mxDouble* y_out = my_mex_get_doubles(y_output_array);
  mxDouble* s_out;
  if (sensitivity || parameter_sensitivity) {
    s_out = my_mex_get_doubles(mex_output[OUTPUT_SENSITIVITY]);
  }
  
  copy_real_to_double(t_values, t_out,  n_points);
  copy_real_to_double(y_data,   y_out, NEQ      );
  
  double* y_out_ptr = y_out + NEQ;



  realtype t          = 0;
  int      rootsfound = false;
  int      i;
  
  for(i = 1; i < n_points; i++) {
    #ifdef PRINT_TIME
    mex_eval("print_temp('t=%.15e')", (double)t_values[i]);
    #endif
    int flag = CVode(cvode, t_values[i], solver_y, &t, CV_NORMAL);
    check(flag, "CVode");
    if (sensitivity || parameter_sensitivity) {
      flag = CVodeGetSens(cvode, &t, solver_s);
      check(flag, "CVodeGetSens");
      copy_real_to_double(s_data, s_out, NEQ * n_sensitivity);
      // todo: fix: s_out seems to be overwritten n_points times here
      // or implement sensitivity ouput at intermediate time values
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
  #ifdef PRINT_TIME
  mex_eval("print_temp()");
  #endif
  
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
    PrintFinalStats(cvode, sensitivity || parameter_sensitivity);
  }
  
  /* Free working memory. Note that matlab would take care of it automatically 
   * if we do not do it here, but freeing memory as soon as we can is better. */
  
  mxDestroyArray(y_output_array);
  mxDestroyArray(t_output_array);

  N_VDestroy(solver_y);
  if (sensitivity || parameter_sensitivity) {
    N_VDestroyVectorArray(solver_s, n_sensitivity);
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
  if (pointer == NULL) {
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

N_Vector get_N_Vector(Argument arg, int size) {
  realtype* data = get_reals(arg, size);
  return N_VMake_Serial(size, data);
}



/*
 * Print some final statistics located in the CVODES memory
 */

void PrintFinalStats(void *cvode, booleantype sensi) {
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
    flag = CVodeGetSensNumErrTestFails(cvode, &netfS);
    check(flag, "CVodeGetSensNumErrTestFails");
    if ((SENSITIVITY_METHOD == CV_STAGGERED) || 
        (SENSITIVITY_METHOD == CV_STAGGERED1)) {
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

realtype* get_reals(Argument arg, int size) {
  check_double(arg);
  check_real  (arg);
  check_size  (arg, size);
  
  realtype* real_data   = malloc(size * sizeof(realtype));
  double*   matlab_data = my_mex_get_doubles(arg.value);
  for (int i = 0; i < size; i++) {
    real_data[i] = (realtype) matlab_data[i];
  }
  return real_data;
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