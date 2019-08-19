#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mex.h>

#include <cvodes/cvodes.h>
#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_band.h>  /* access to band SUNMatrix       */
#include <sunlinsol/sunlinsol_band.h>  /* access to band SUNLinearSolver */
#include <sunmatrix/sunmatrix_dense.h>  /* access to dense SUNMatrix       */
#include <sunlinsol/sunlinsol_dense.h>  /* access to dense SUNLinearSolver */

#include <sundials/sundials_types.h>
#include <sundials/sundials_math.h>

#include "user_data.h"




#define output_t           0
#define output_y           1
#define output_sensitivity 2


/* Problem Constants */


#define ABS_TOL RCONST(1.e-6) /* default scalar absolute tolerance */
#define REL_TOL RCONST(1.e-6) /* default relative tolerance */



#define NS    1

#define INITIAL_OUTPUT_SIZE

#define INCREASING 1

#define ZERO  RCONST(0.0)

#define MAX_NUM_STEPS 10*1000*1000

#define debug false

/* Functions Called by the CVODES Solver */


// implemented in dydt_cvode.c
int dydt_ode(realtype t, N_Vector u, N_Vector udot, void *);

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

static int return_to_plane(realtype t, N_Vector u, realtype* g  , void *);


/* Private Helper Functions */


static void check_null(void* ptr, const char* funcname);
static void check_return_value(int return_value, const char* funcname);

static realtype* my_mex_get_doubles(const mxArray* array);
static void check_double(const mxArray* array,               char* arrayname);
static void check_bool  (const mxArray* array,               char* arrayname);
static void check_size  (const mxArray* array, int expected, char* arrayname);
static N_Vector get_N_Vector(const mxArray* array, char* input_name, int size);
static void PrintFinalStats(void *cvode_mem, booleantype sensi,
                            booleantype err_con, int sensi_meth);

static int         nPoints;
static realtype*   t_values;
static N_Vector    initial_point;
static realtype*   parameters;
static int         max_num_points;
static booleantype max_num_points_given;

static booleantype sensitivity;
static N_Vector    sensitivity_vector;

static booleantype cycle_detection;
static realtype    lower_bound_period;
static booleantype lower_bound_period_given;
static N_Vector    point_on_limitcycle;
static N_Vector    tangent_to_limitcycle;

static realtype    abs_tol;
static booleantype abs_tol_given;

static realtype    rel_tol;
static booleantype rel_tol_given;

void mexFunction(int n_output,       mxArray *mex_output[], 
                 int n_input,  const mxArray *mex_input []  ) {
  
  nPoints                  = 0;
  t_values                 = NULL;
  initial_point            = NULL;
  parameters               = NULL;
  max_num_points           = 0;
  max_num_points_given     = false;

  sensitivity              = false;
  sensitivity_vector       = NULL;

  cycle_detection          = false;
  lower_bound_period       = 0;
  lower_bound_period_given = false;
  point_on_limitcycle      = NULL;
  tangent_to_limitcycle    = NULL;
  
  abs_tol                  = ABS_TOL;
  abs_tol_given            = false;
  
  rel_tol                  = REL_TOL;
  rel_tol_given            = false;
  
  UserData user_data = malloc(sizeof(UserData));
  check_null(user_data, "malloc user_data");

  

  
  /* check for proper number of inputs */
  if(n_input < 6) {
    mexErrMsgIdAndTxt("cvode:n_input", 
                         "at least 4 name - value pairs are required as input");
  }
  
  if (n_input % 2 != 0) {
    mexErrMsgIdAndTxt("cvode:nrhs_uneven",
            "Error: the number of arguments is uneven. "
            "The number for arguments is %d. "
            "Arguments must be supplied as name-value pairs", n_input);
  }
    
  
  for ( int i = 0; i < n_input; i += 2 ) {
    
    
    if ( !mxIsChar(mex_input[i]) ) {
      mexErrMsgIdAndTxt("cvode:not_name_value", 
                              "arguments must be supplied in name-value pairs");
    }
    char* arg_name = (char*) mxArrayToString(mex_input[i]);

    if        ( ! strcmp(arg_name, "t_values"            ) ) {
      
      check_double(mex_input[i + 1], "t_values");
      t_values = (realtype*) my_mex_get_doubles(mex_input[i + 1]);
      
      nPoints      = mxGetNumberOfElements(mex_input[i + 1]);
      if (nPoints < 2) {
        mexErrMsgIdAndTxt("cvode:wrong_size",
                        "Input vector t_values must have at least 2 elements.");
      }
      
    } else if ( ! strcmp(arg_name, "initial_point"       ) ) {
      
      initial_point = get_N_Vector(mex_input[i + 1], "initial_point", NEQ);
      
    } else if ( ! strcmp(arg_name, "ode_parameters"       ) ) {
      
      check_double(mex_input[i + 1],               "ode_parameters");
      check_size  (mex_input[i + 1], n_parameters, "ode_parameters");
      
      parameters = my_mex_get_doubles(mex_input[i + 1]);
      user_data->parameters = my_mex_get_doubles(mex_input[i + 1]);
      
    } else if ( ! strcmp(arg_name, "sensitivity_vector"   ) ) {
      
      sensitivity_vector
                    = get_N_Vector(mex_input[i + 1], "sensitivity_vector", NEQ);
      sensitivity = true;
      
    } else if ( ! strcmp(arg_name, "cycle_detection"      ) ) {
      
      check_bool(mex_input[i + 1], "cycle_detection");
      
      cycle_detection = mxIsLogicalScalarTrue(mex_input[i + 1]);
      
    } else if ( ! strcmp(arg_name, "lower_bound_period"  ) ) {
      
      check_double(mex_input[i + 1],    "lower_bound_period");
      check_size  (mex_input[i + 1], 1, "lower_bound_period");
      
      lower_bound_period = mxGetScalar(mex_input[i + 1]);
      lower_bound_period_given = true;
      
    } else if ( ! strcmp(arg_name, "point_on_limitcycle"  ) ) {
      
      point_on_limitcycle = 
                     get_N_Vector(mex_input[i + 1], "point_on_limitcycle", NEQ);
      tangent_to_limitcycle = N_VNew_Serial(NEQ);
      
    } else if ( ! strcmp(arg_name, "abs_tol"  ) ) {
      
      check_double(mex_input[i + 1],    "abs_tol");
      check_size  (mex_input[i + 1], 1, "abs_tol");
      abs_tol = mxGetScalar(mex_input[i + 1]);
      if (abs_tol <= 0) {
        mexErrMsgIdAndTxt("cvode:abs_tol_not_positive", 
                                         "Error: abs_tol is not positive");   
      }
      abs_tol_given = true;
      
    } else if ( ! strcmp(arg_name, "rel_tol"  ) ) {
      
      check_double(mex_input[i + 1],    "rel_tol");
      check_size  (mex_input[i + 1], 1, "rel_tol");
      abs_tol = mxGetScalar(mex_input[i + 1]);
      if (abs_tol <= 0) {
        mexErrMsgIdAndTxt("cvode:abs_tol_not_positive", 
                                         "Error: rel_tol is not positive");   
      }
      rel_tol_given = true;
      
    } else {
  
      mexErrMsgIdAndTxt("cvode:invalid_argument_name", 
                        "Error: Invalid argument name: %s.", arg_name);
    }
  }
  
  
  if(n_output <  2) {
    mexErrMsgIdAndTxt("cvode:n_output", "2 outputs are needed: t and y "
                     "(or 3 outputs for sensitivity analysis)");
  }
  
  if(sensitivity && n_output < 3) {
    mexErrMsgIdAndTxt("cvode:n_output",
                      "Error: You enabled sensitivity analysis, " 
                      "but you only specified less than three outputs. "
                      "Three outputs are needed for sensitivity analysis: "
                      "t, y, and s");
  }

  if (t_values == NULL) {
    mexErrMsgIdAndTxt("cvode:no_t_values", 
                                         "Error: you did not specify t_values");
  }
  
  if (initial_point == NULL) {
    mexErrMsgIdAndTxt("cvode:no_y_values",
                                    "Error: you did not specify initial_point");
  }
  
  if (parameters == NULL) {
    mexErrMsgIdAndTxt("cvode:no_params",
                                   "Error: you did not specify ode_parameters");
  }

  if (cycle_detection) {
    // we wait to compute tangent_to_limitcycle, since we need parameters to 
    // be loaded first
    dydt_ode(0, point_on_limitcycle, tangent_to_limitcycle, user_data);
    
    if (sensitivity) {
       mexErrMsgIdAndTxt("cvode:sens_and_cd",
                             "Error: simultaneous cycle detection and "
                             "sensitivity analysis is not supported");   
    }
    if (! lower_bound_period_given) {
      mexErrMsgIdAndTxt("cvode:no_lower_bound_period",
                             "Error: You enabled cycle detection, "
                             "but you did not specify lower_bound_period");
    }
    if ( lower_bound_period <= 0 ) {
      mexErrMsgIdAndTxt("cvode:gap_tolerance_not_positive",
                             "Error: You enabled cycle detection, "
                             "but lower_bound_period is not positive");
    }
    if (point_on_limitcycle == NULL) {
      mexErrMsgIdAndTxt("cvode:no_point_on_limitcycle",
                             "Error: You enabled cycle detection, "
                             "but you did not specify point_on_limit_cycle");    
    }
  }
  
  if (sensitivity) {
    if (n_output != 3) {
      mexErrMsgIdAndTxt("cvode:n_output",
                             "Error: you enabled sensitivity analysis "
                             "but you did not specify three outputs");   
    }
  }

  void*             cvode_mem   = NULL;
  SUNMatrix         A           = NULL;
  SUNLinearSolver   LS          = NULL;
  int               sensi_meth  = CV_SIMULTANEOUS;
  booleantype       err_con     = false;

  
  N_Vector y = N_VNew_Serial(NEQ);  
  N_Vector s = N_VNew_Serial(NEQ);
  
  double* initial_point_data = N_VGetArrayPointer(initial_point);
  double* y_data             = N_VGetArrayPointer(y);
  memcpy(y_data, initial_point_data, NEQ * sizeof(double));
  
  
  double* s_data;
  if (sensitivity) {
    double* sens_vector_data = N_VGetArrayPointer(sensitivity_vector);
    s_data                   = N_VGetArrayPointer(s); 
    memcpy(s_data, sens_vector_data, NEQ * sizeof(double));
  }
  
  
  
  

  
  /* Create CVODES object */
  cvode_mem = CVodeCreate(CV_BDF);
  check_null((void *)cvode_mem, "CVodeCreate");

  /* Allocate CVODES memory */
  int retval;
  retval = CVodeInit(cvode_mem, dydt_ode, t_values[0], y);
  check_return_value(retval, "CVodeInit");
  
  retval = CVodeSetUserData(cvode_mem, user_data);
  check_return_value(retval, "CVodeSetUserData");

  retval = CVodeSStolerances(cvode_mem, rel_tol, abs_tol);
  check_return_value(retval, "CVodeSStolerances");
  

  
  /* Create banded SUNMatrix for use in linear solves 
   * -- since this will be factored, 
   *  set the storage bandwidth to be the sum of upper and lower bandwidths */
  
  #if ANALYTIC_JACOBIAN
  A = SUNDenseMatrix(NEQ,NEQ);
  check_null((void *)A, "SUNDenseMatrix");
  #else
  A = SUNBandMatrix(NEQ, 10, 10);
  check_null((void *)A, "SUNBandMatrix");
  #endif
  /* Create banded SUNLinearSolver object for use by CVode */
  
  #if ANALYTIC_JACOBIAN
  LS = SUNLinSol_Dense(y,A);//SUNLinSol_Band(u, A);
  check_null((void *)LS, "SUNLinSol_Dense") ;
  #else
  LS = SUNLinSol_Band(y, A);
  check_null((void *)LS, "SUNLinSol_Band") ;
  #endif
    /* Call CVodeSetLinearSolver to attach the matrix and linear solver to CVode */
  retval = CVodeSetLinearSolver(cvode_mem, LS, A);
  check_return_value(retval, "CVodeSetLinearSolver");
  
    /* Set the user-supplied Jacobian routine Jac */
  #if ANALYTIC_JACOBIAN
  retval = CVodeSetJacFn(cvode_mem, jacobian_dydt);
  #endif
  
  check_return_value(retval, "CVodeSetJacFn");
  
  CVodeSetMaxNumSteps(cvode_mem, MAX_NUM_STEPS);
  
  int rootfinding_direction = INCREASING;
  if (cycle_detection) {
    CVodeRootInit(cvode_mem, 1, return_to_plane);
    CVodeSetRootDirection(cvode_mem, &rootfinding_direction);
  }
 

  
  if(sensitivity) {

    retval = CVodeSensInit1(cvode_mem, NS, sensi_meth, d_sensitivity_dt, &s);
    //retval = CVodeSensInit1(cvode_mem, NS, sensi_meth, NULL, &uS);
    check_return_value(retval, "CVodeSensInit1");
    CVodeSensSStolerances(cvode_mem, rel_tol, &abs_tol);
    //retval = CVodeSensEEtolerances(cvode_mem);
    //check_return_value(retval, "CVodeSensEEtolerances");

    retval = CVodeSetSensErrCon(cvode_mem, err_con);
    check_return_value(retval, "CVodeSetSensErrCon");

    retval = CVodeSetSensDQMethod(cvode_mem, CV_CENTERED, ZERO);
    check_return_value(retval, "CVodeSetSensDQMethod");
    
    realtype dummy_value = 0;
    retval = CVodeSetSensParams(cvode_mem, &dummy_value, NULL, NULL);
    check_return_value(retval, "CVodeSetSensParams");
    if (debug) {
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
    if (debug) mexPrintf("Sensitivity: NO \n");
  }
  
  mex_output[output_t]     = mxCreateDoubleMatrix(1  , nPoints, mxREAL);
  mxArray* y_output_buffer = mxCreateDoubleMatrix(NEQ, nPoints, mxREAL);
  if (sensitivity) {
    mex_output[output_sensitivity] = mxCreateDoubleMatrix(NEQ, 1, mxREAL);
  }
  
  mxDouble* t_out = my_mex_get_doubles(mex_output[output_t]);
  mxDouble* y_out = my_mex_get_doubles(y_output_buffer);
  mxDouble* s_out;
  if (sensitivity) {
    s_out = my_mex_get_doubles(mex_output[output_sensitivity]);
  }
  
  memcpy(t_out, t_values, nPoints * sizeof(double));
  memcpy(y_out, y_data,   NEQ     * sizeof(double));

  realtype* y_out_ptr = y_out + NEQ;



  realtype t=0;
  int rootsfound = false;
 
  int i;
  
  for(i = 1; i < nPoints; i++) {  
    retval = CVode(cvode_mem, t_values[i], y, &t, CV_NORMAL);
    check_return_value(retval, "CVode");
    if (sensitivity) {
      retval = CVodeGetSens(cvode_mem, &t, &s);
      check_return_value(retval, "CVodeGetSens");
      
      s_data = N_VGetArrayPointer(s);
      memcpy(s_out, s_data, NEQ * sizeof(double));  
    }
    if (cycle_detection) {
      CVodeGetRootInfo(cvode_mem, &rootsfound);
    }
    
    memcpy(y_out_ptr, y_data, NEQ * sizeof(double));
    y_out_ptr += NEQ;

    if (rootsfound) {
      //mexPrintf("cvode: period: %.6f\n", t);
      t_out[i] = t;
      break;
    }
  }
  if (i+1 < nPoints) {
    mxSetN(mex_output[output_t], i+1);
    mxSetN(y_output_buffer, i+1);
  }
  mxArray* transposed = mxCreateDoubleMatrix(i+1, NEQ, mxREAL);
  mexCallMATLAB(1, &transposed, 1, &y_output_buffer, "transpose");
  mex_output[output_y] = transposed;
   
   PrintFinalStats(cvode_mem, sensitivity, err_con, sensi_meth);
  
  /* Free memory */
  
  mxDestroyArray(y_output_buffer);

  N_VDestroy(y);
  N_VDestroy(s);

  CVodeFree(&cvode_mem);


}

static void check_null(void* ptr, const char* funcname) {
  if (ptr == NULL) {
    mexErrMsgIdAndTxt("CVODE:null_pointer",
             "SUNDIALS_ERROR: %s() failed - returned NULL pointer\n", funcname); 
  }
}

static void check_return_value(int return_value, const char* funcname) {
  if (return_value < 0) {
    mexErrMsgIdAndTxt("CVODE:error", 
                      "SUNDIALS_ERROR: %s() failed with error_code %s",
                      funcname, CVodeGetReturnFlagName(return_value));
  }
}

realtype* my_mex_get_doubles(const mxArray* array) {
  #if MX_HAS_INTERLEAVED_COMPLEX
  return (realtype*) mxGetDoubles(array);
  #else
  return (realtype*) mxGetPr(     array);
  #endif 
}

static int return_to_plane(realtype t, N_Vector y, realtype *gout, 
                                                              void *user_data) {
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

static void check_double(const mxArray* array, char* arrayname) {
  if ( !mxIsDouble(array) ) {
    mexErrMsgIdAndTxt("cvode:not_double", 
               "Input vector %s is not a vector of doubles.", arrayname);
  }
}

static void check_bool(const mxArray* array, char* arrayname) {
  if ( !mxIsLogicalScalar(array) ) {
    mexErrMsgIdAndTxt("cvode:not_bool", 
                                       "Input %s is not a boolean.", arrayname);
  }
}


static void check_size(const mxArray* array, int size, char* arrayname) {
  if ( mxGetNumberOfElements(array) != size ) {
    mexErrMsgIdAndTxt("cvode:incorrect_size",  
            "Error: Input vector %s does not have the correct size. "
            "The correct size is %d. The actual size is %d.",
            arrayname, size, (int)mxGetNumberOfElements(array));
  }
}

static N_Vector get_N_Vector(const mxArray* array, char* input_name, int size) {
  check_double(array,       input_name);
  check_size  (array, size, input_name);
  realtype* data = my_mex_get_doubles(array);
  return N_VMake_Serial(size, data);
}


/*
 * Print some final statistics located in the CVODES memory
 */

static void PrintFinalStats(void *cvode_mem, booleantype sensi,
                            booleantype err_con, int sensi_meth)
{
  long int nst;
  long int nfe, nsetups, nni, ncfn, netf;
  long int nfSe, nfeS, nsetupsS, nniS, ncfnS, netfS;
  int retval;

  retval = CVodeGetNumSteps(cvode_mem, &nst);
  //check_retval(&retval, "CVodeGetNumSteps", 1);
  retval = CVodeGetNumRhsEvals(cvode_mem, &nfe);
  //check_retval(&retval, "CVodeGetNumRhsEvals", 1);
  retval = CVodeGetNumLinSolvSetups(cvode_mem, &nsetups);
  //check_retval(&retval, "CVodeGetNumLinSolvSetups", 1);
  retval = CVodeGetNumErrTestFails(cvode_mem, &netf);
  //check_retval(&retval, "CVodeGetNumErrTestFails", 1);
  retval = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
  //check_retval(&retval, "CVodeGetNumNonlinSolvIters", 1);
  retval = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
  //check_retval(&retval, "CVodeGetNumNonlinSolvConvFails", 1);

  if (sensi) {
    retval = CVodeGetSensNumRhsEvals(cvode_mem, &nfSe);
    //check_retval(&retval, "CVodeGetSensNumRhsEvals", 1);
    retval = CVodeGetNumRhsEvalsSens(cvode_mem, &nfeS);
    //check_retval(&retval, "CVodeGetNumRhsEvalsSens", 1);
    retval = CVodeGetSensNumLinSolvSetups(cvode_mem, &nsetupsS);
    //check_retval(&retval, "CVodeGetSensNumLinSolvSetups", 1);
    if (err_con) {
      retval = CVodeGetSensNumErrTestFails(cvode_mem, &netfS);
      //check_retval(&retval, "CVodeGetSensNumErrTestFails", 1);
    } else {
      netfS = 0;
    }
    if ((sensi_meth == CV_STAGGERED) || (sensi_meth == CV_STAGGERED1)) {
      retval = CVodeGetSensNumNonlinSolvIters(cvode_mem, &nniS);
      //check_retval(&retval, "CVodeGetSensNumNonlinSolvIters", 1);
      retval = CVodeGetSensNumNonlinSolvConvFails(cvode_mem, &ncfnS);
      //check_retval(&retval, "CVodeGetSensNumNonlinSolvConvFails", 1);
    } else {
      nniS = 0;
      ncfnS = 0;
    }
  }

  mexPrintf("\nFinal Statistics\n\n");
  mexPrintf("number of steps     = %5ld\n\n", nst);
  mexPrintf("number of rhs evals     = %5ld\n",   nfe);
  mexPrintf("netf    = %5ld    nsetups  = %5ld\n", netf, nsetups);
  mexPrintf("nni     = %5ld    ncfn     = %5ld\n", nni, ncfn);

  if(sensi) {
    mexPrintf("\n");
    mexPrintf("nfSe    = %5ld    nfeS     = %5ld\n", nfSe, nfeS);
    mexPrintf("netfs   = %5ld    nsetupsS = %5ld\n", netfS, nsetupsS);
    mexPrintf("nniS    = %5ld    ncfnS    = %5ld\n", nniS, ncfnS);
  }

}