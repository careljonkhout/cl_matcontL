#include <stdio.h>
#include <stdlib.h>
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

#define OUTPUT_T           0
#define OUTPUT_Y           1
#define OUTPUT_SENSITIVITY 2
#define DEFAULT_ABS_TOL    RCONST(1.e-6) /* default scalar absolute tolerance */
#define DEFAULT_REL_TOL    RCONST(1.e-6) /* default relative tolerance */
#define NS                 1
#define INCREASING         1
#define MAX_NUM_STEPS      10*1000*1000

// implemented in dydt_cvode.c
int dydt_cvode(double t, N_Vector u, N_Vector udot, void *);

#if ANALYTIC_JACOBIAN
// implemented in jacobian_cvode.c
int jacobian_dydt(
               double t, N_Vector y_vec, N_Vector fy, SUNMatrix jac_structure,
              void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
#endif
// implemented in d_sensitivity_dt.c
int d_sensitivity_dt(int Ns, double t,
                N_Vector y_vector, 
                N_Vector ydot,
                int iS,
                N_Vector s_vector, 
                N_Vector ds_vector,
                void *user_data,
                N_Vector tmp1, N_Vector tmp2);


int return_to_plane(double t, N_Vector u, double* g  , void *);

void error_handler(int error_code, const char *module,
                          const char *function, char *msg, void *eh_data);


double* my_mex_get_doubles(const mxArray* array);
void check_double    (const mxArray* array,               char* arrayname);
void check_bool      (const mxArray* array,               char* arrayname);
void check_size      (const mxArray* array, int expected, char* arrayname);
void check_positive  (const mxArray* array,               char* arrayname);
N_Vector get_N_Vector(const mxArray* array, char* input_name, int size);
double*  get_doubles (const mxArray* array, char* input_name, int size);
N_Vector new_N_Vector(int length);

void check_argument_given(void* pointer,                char* arrayname);

void PrintFinalStats(void *cvode, booleantype sensi,
                            booleantype err_con, int sensi_meth);


double       lower_bound_period;
N_Vector     point_on_limitcycle;
N_Vector     tangent_to_limitcycle;

void mexFunction(int n_output,       mxArray *mex_output[], 
                 int n_input,  const mxArray *mex_input []  ) {
  
  int          n_points                 = 0;
  double*      t_values                 = NULL;
  double*      initial_point            = NULL;

  int          sensitivity              = false;
  double*      sensitivity_vector       = NULL;

  booleantype  cycle_detection          = false;
               lower_bound_period       = 0;
  booleantype  lower_bound_period_given = false;
               point_on_limitcycle      = NULL;
               tangent_to_limitcycle    = NULL; 

  booleantype  verbose                  = false;
  
  double       abs_tol                  = DEFAULT_ABS_TOL;
  double       rel_tol                  = DEFAULT_REL_TOL;
  /* abs_tol and rel_tol will be overwritten if the caller specifies them */
  
  
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
    char* arg_name = (char*) mxArrayToString(mex_input[i]);

    if        ( ! strcmp(arg_name, "t_values"             ) ) {
      
      check_double(mex_input[i + 1], (char*) "t_values");
      t_values = (double*) my_mex_get_doubles(mex_input[i + 1]);
      
      n_points      = mxGetNumberOfElements(mex_input[i + 1]);
      if (n_points < 2) {
        mexErrMsgIdAndTxt("cvode:wrong_size",
                        "Input vector t_values must have at least 2 elements.");
      }
      
    } else if ( ! strcmp(arg_name, "initial_point"        ) ) {
      
      initial_point = get_doubles(mex_input[i + 1], "initial_point", NEQ);
      
    } else if ( ! strcmp(arg_name, "ode_parameters"       ) ) {
      
      check_double(mex_input[i + 1],               "ode_parameters");
      check_size  (mex_input[i + 1], N_PARAMETERS, "ode_parameters");
      
      user_data->parameters = my_mex_get_doubles(mex_input[i + 1]);
      
    } else if ( ! strcmp(arg_name, "sensitivity_vector"   ) ) {
      
      sensitivity_vector
                     = get_doubles(mex_input[i + 1], "sensitivity_vector", NEQ);
      
      sensitivity = true;
      
    } else if ( ! strcmp(arg_name, "cycle_detection"      ) ) {
      
      check_bool(mex_input[i + 1], "cycle_detection");
      
      cycle_detection = mxIsLogicalScalarTrue(mex_input[i + 1]);
      
    } else if ( ! strcmp(arg_name, "lower_bound_period"  ) ) {
      
      check_double  (mex_input[i + 1],    "lower_bound_period");
      check_size    (mex_input[i + 1], 1, "lower_bound_period");
      check_positive(mex_input[i + 1],    "lower_bound_period");
      
      lower_bound_period = mxGetScalar(mex_input[i + 1]);
      lower_bound_period_given = true;
      
    } else if ( ! strcmp(arg_name, "point_on_limitcycle"  ) ) {
      
      point_on_limitcycle = 
                     get_N_Vector(mex_input[i + 1], "point_on_limitcycle", NEQ);
      tangent_to_limitcycle = new_N_Vector(NEQ);
      
    } else if ( ! strcmp(arg_name, "abs_tol"  ) ) {
      
      check_double  (mex_input[i + 1],    "abs_tol");
      check_size    (mex_input[i + 1], 1, "abs_tol");
      check_positive(mex_input[i + 1],    "abs_tol");
      abs_tol = mxGetScalar(mex_input[i + 1]);
      
    } else if ( ! strcmp(arg_name, "rel_tol"  ) ) {
      
      check_double  (mex_input[i + 1],    "rel_tol");
      check_size    (mex_input[i + 1], 1, "rel_tol");
      check_positive(mex_input[i + 1],    "rel_tol");
      rel_tol = mxGetScalar(mex_input[i + 1]);
      
    } else if ( ! strcmp(arg_name, "verbose"      ) ) {
      
      check_bool(mex_input[i + 1], "verbose");
      verbose = mxIsLogicalScalarTrue(mex_input[i + 1]);
      
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

  
  double* y_data             = N_VGetArrayPointer(solver_y);
  
  memcpy(y_data, initial_point, NEQ * sizeof(double));
  
  double* s_data;
  if ( sensitivity ) {
    s_data = N_VGetArrayPointer(solver_s); 
    memcpy(s_data, sensitivity_vector, NEQ * sizeof(double));
  }
  
  cvode = CVodeCreate_nullchecked(CV_BDF);

  CVodeInit_checked           (cvode, dydt_cvode    , t_values[0], solver_y);  
  CVodeSetErrHandlerFn_checked(cvode, error_handler , NULL                 );
  CVodeSetUserData_checked    (cvode, user_data                            );
  CVodeSStolerances_checked   (cvode, rel_tol       , abs_tol              );
  CVodeSetMaxNumSteps_checked (cvode, MAX_NUM_STEPS                        );
     
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
    double dummy_value = 0;
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
  
  memcpy(t_out, t_values, n_points * sizeof(double));
  memcpy(y_out, y_data,   NEQ      * sizeof(double));

  double* y_out_ptr = y_out + NEQ;



  double t=0;
  int rootsfound = false;
 
  int i;
  
  int n_points_to_compute = sensitivity ? 2 : n_points;
  
  for(i = 1; i < n_points; i++) {  
    int flag = CVode(cvode, t_values[i], solver_y, &t, CV_NORMAL);
    check(flag, "CVode");
    if (sensitivity) {
      flag = CVodeGetSens(cvode, &t, &solver_s);
      check(flag, "CVodeGetSens");
      memcpy(s_out, s_data, NEQ * sizeof(double));  
    }
    if (cycle_detection) {
      CVodeGetRootInfo(cvode, &rootsfound);
    }
    
    memcpy(y_out_ptr, y_data, NEQ * sizeof(double));
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


double* my_mex_get_doubles(const mxArray* array) {
  #if MX_HAS_INTERLEAVED_COMPLEX
  return (double*) mxGetDoubles(array);
  #else
  return (double*) mxGetPr(     array);
  #endif 
}

int return_to_plane(double t, N_Vector y, double *gout, 
                                                              void *user_data) {
  if (t > lower_bound_period) {
    *gout = 0;
    double* y_data   = NV_DATA_S(y);
    double* pol_data = NV_DATA_S(point_on_limitcycle);
    double* ttl_data = NV_DATA_S(tangent_to_limitcycle);
    for (int i = 0; i < NV_LENGTH_S(y); i++) {
      *gout += (y_data[i] - pol_data[i]) * ttl_data[i];
    }
  } else {
    *gout = 1;
  }
  return 0;
}

#define BUFFER_LENGTH 1000

void error_handler(int error_code, const char* module,
         const char* function, char* message, void *eh_data) {
  if (error_code < 0) {
    mexErrMsgIdAndTxt("cvode:integrator_error",
            "Error in %s() in sundails module %s. "
            "error code: %d: %s.\n  %s\n",
            function, module,
            error_code, CVodeGetReturnFlagName(error_code), message);
  } else {
    char* message_buffer = mxMalloc(BUFFER_LENGTH * sizeof(char));
    char* buffer = mxMalloc((BUFFER_LENGTH - 100) * sizeof(char));
    snprintf(message_buffer, BUFFER_LENGTH - 100, 
            "Error in %s() in sundails module %s. error code: %d: %s. %s", 
            function, module, error_code, 
            CVodeGetReturnFlagName(error_code), message);       
    snprintf(buffer, BUFFER_LENGTH, "print_diag(0,'%s\\\n')", message_buffer);
    mexEvalString(buffer);
  }
}

void check_double(const mxArray* array, char* arrayname) {
  if ( !mxIsDouble(array) ) {
    mexErrMsgIdAndTxt("cvode:not_double", 
             "Input vector %s is not a vector of doubles.", arrayname);
  }
}

void check_bool(const mxArray* array, char* arrayname) {
  if ( !mxIsLogicalScalar(array) ) {
    mexErrMsgIdAndTxt("cvode:not_bool", 
            "Input %s is not a boolean.", arrayname);
  }
}

void check_size(const mxArray* array, int size, char* arrayname) {
  if ( mxGetNumberOfElements(array) != size ) {
    mexErrMsgIdAndTxt("cvode:incorrect_size",  
            "Error: Input vector %s does not have the correct size. "
            "The correct size is %d. The actual size is %d.",
            arrayname, size, (int)mxGetNumberOfElements(array));
  }
}

void check_positive(const mxArray* scalar, char* scalarname) {
  if (mxGetScalar(scalar) <= 0) {
    mexErrMsgIdAndTxt("cvode:not_positive",  
            "Error: Input %s is not positive.", scalarname);
  }
}

void check_argument_given(void * pointer, char* argument_name) {
  if (argument_name == NULL) {
    mexErrMsgIdAndTxt("cvode:missing_argument",  
                      "Argument '%s' was not passed to cvode.", argument_name);  
  }
}

double* get_doubles(const mxArray* array, char* input_name, int size) {
  check_double(array,       input_name);
  check_size  (array, size, input_name);
  return my_mex_get_doubles(array);
}

N_Vector get_N_Vector(const mxArray* array, char* input_name, int size) {
  check_double(array,       input_name);
  check_size  (array, size, input_name);
  double* data = my_mex_get_doubles(array);
  return N_VMake_Serial(size, data);
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