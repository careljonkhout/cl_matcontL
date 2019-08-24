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
#include "cvode_with_checks.c"



#define output_t           0
#define output_y           1
#define output_sensitivity 2
#define ABS_TOL            RCONST(1.e-6) /* default scalar absolute tolerance */
#define REL_TOL            RCONST(1.e-6) /* default relative tolerance */
#define NS                 1
#define INCREASING         1
#define ZERO               RCONST(0.0)
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

static int dydt_simultaneous(
        realtype t, N_Vector u_and_s, N_Vector dydt_vec, void *user_data);

static int jacobian_simultaneous(realtype t, N_Vector u_and_s, N_Vector fy,
              SUNMatrix jacobian, void* user_data,
              N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

static int return_to_plane(realtype t, N_Vector u, realtype* g  , void *);

static void error_handler(int error_code, const char *module,
                          const char *function, char *msg, void *eh_data);

static realtype* my_mex_get_doubles(const mxArray* array);
static void check_double(const mxArray* array,               char* arrayname);
static void check_bool  (const mxArray* array,               char* arrayname);
static void check_size  (const mxArray* array, int expected, char* arrayname);
static N_Vector get_N_Vector(const mxArray* array, char* input_name, int size);
// static void append_to_array(mxArray* array, int from_index,
//                               double* data, int data_length);
// static void grow_mxArray(mxArray* array);
static void PrintFinalStats(void *cvode, booleantype sensi,
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

static booleantype verbose;

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
  
  verbose                  = false;
  
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
      check_size  (mex_input[i + 1], N_PARAMETERS, "ode_parameters");
      
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
      
    } else if ( ! strcmp(arg_name, "verbose"      ) ) {
      
      check_bool(mex_input[i + 1], "verbose");
      
      verbose = mxIsLogicalScalarTrue(mex_input[i + 1]);
      
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
    dydt_cvode(0, point_on_limitcycle, tangent_to_limitcycle, user_data);
    
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

  void*             cvode   = NULL;
  SUNMatrix         A           = NULL;
  SUNLinearSolver   LS          = NULL;
  int               sensi_meth  = CV_SIMULTANEOUS;
  booleantype       err_con     = true;

  int size = (MANUAL_SENSITIVITY && sensitivity) ? 2 * NEQ : NEQ;
  
  N_Vector y = N_VNew_Serial(size); 
  
  N_Vector s;
  if (sensitivity && !MANUAL_SENSITIVITY) {
    s = N_VNew_Serial(NEQ);
  }
  
  double* initial_point_data = N_VGetArrayPointer(initial_point);
  double* y_data             = N_VGetArrayPointer(y);
  memcpy(y_data, initial_point_data, NEQ * sizeof(double));
  
  
  double* s_data;
  if (sensitivity) {
    double* sens_vector_data = N_VGetArrayPointer(sensitivity_vector);
    if (MANUAL_SENSITIVITY) {
      memcpy(y_data + NEQ, sens_vector_data, NEQ * sizeof(double));
    } else {
      s_data                   = N_VGetArrayPointer(s); 
      memcpy(s_data      , sens_vector_data, NEQ * sizeof(double));
    }
  }
  
  cvode = CVodeCreate_nullchecked(CV_BDF);

  if (sensitivity && MANUAL_SENSITIVITY) {
    CVodeInit_checked(cvode, dydt_simultaneous, t_values[0], y);  
  } else {
    CVodeInit_checked(cvode, dydt_cvode       , t_values[0], y);
  }
  CVodeSetErrHandlerFn_checked(cvode, error_handler , NULL          );
  CVodeSetUserData_checked    (cvode, user_data                     );
  CVodeSStolerances_checked   (cvode, rel_tol       , abs_tol       );
  CVodeSetMaxNumSteps_checked (cvode, MAX_NUM_STEPS                 );
     
  #if DENSE_JACOBIAN
  A  = SUNDenseMatrix_nullchecked(size, size);
  LS = SUNLinSol_Dense_nullchecked(y, A);
  #endif
  
  #if BANDED_JACOBIAN
  A  = SUNBandMatrix_nullchecked(size, 2, 2);
  LS = SUNLinSol_Band_nullchecked(y, A);
  #endif
   
  CVodeSetLinearSolver_checked(cvode, LS, A);
 
  /* Set the user-supplied Jacobian routine Jac */
  #if ANALYTIC_JACOBIAN
  if (sensitivity && MANUAL_SENSITIVITY) {
    CVodeSetJacFn(cvode, jacobian_simultaneous);
  } else {
    CVodeSetJacFn(cvode, jacobian_dydt);
  }
  #endif
  
  if (sensitivity) {
    //CVodeSetMaxStep(cvode, 0.0005);
  }
  
  
  
  int rootfinding_direction = INCREASING;
  if (cycle_detection) {
    CVodeRootInit        (cvode, 1,                     return_to_plane);
    CVodeSetRootDirection(cvode, &rootfinding_direction                );
  }
 

  
  if(sensitivity && ! MANUAL_SENSITIVITY) {
    CVodeSensInit1       (cvode, NS, sensi_meth, d_sensitivity_dt, &s);
    CVodeSensEEtolerances(cvode);
    //CVodeSensSStolerances(cvode, rel_tol, &abs_tol);
    CVodeSetSensErrCon   (cvode, err_con);    
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
  
  
  
  mex_output[output_t]     = mxCreateDoubleMatrix(1  , nPoints, mxREAL);
  
  mxArray* y_output_buffer = mxCreateDoubleMatrix(size, nPoints, mxREAL);
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
    int flag = CVode(cvode, t_values[i], y, &t, CV_NORMAL);
    check(flag, "CVode");
    if (sensitivity && ! MANUAL_SENSITIVITY) {
      flag = CVodeGetSens(cvode, &t, &s);
      check(flag, "CVodeGetSens");
      
      s_data = N_VGetArrayPointer(s);
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
  
  if (sensitivity && MANUAL_SENSITIVITY) {
    memcpy(s_out, y_data + NEQ, NEQ * sizeof(double));  
  }
  
  if (i+1 < nPoints) {
    mxSetN(mex_output[output_t], i+1);
    mxSetN(y_output_buffer     , i+1);
  }
  mxArray* transposed = mxCreateDoubleMatrix(i+1, NEQ, mxREAL);
  mexCallMATLAB(1, &transposed, 1, &y_output_buffer, "transpose");
  mex_output[output_y] = transposed;
   
  if (verbose) {
    PrintFinalStats(cvode, sensitivity, err_con, sensi_meth);
  }
  
  /* Free memory */
  
  mxDestroyArray(y_output_buffer);

  N_VDestroy(y);
  if (sensitivity && !MANUAL_SENSITIVITY) {
    N_VDestroy(s);
  }

  CVodeFree(&cvode);
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

#define X(i) u[2*(i)]
#define Y(i) u[2*(i)+1]

#define SX(i) s[2*(i)]
#define SY(i) s[2*(i)+1]

#define DSX_DT(i) dsdt[2*(i)]
#define DSY_DT(i) dsdt[2*(i)+1]

#define L  parameters[0]
#define A  parameters[1]
#define B  parameters[2]
#define DX parameters[3]
#define DY parameters[4]

void dydt_ode(double* y, double* dydt, double* parameters);

static int dydt_simultaneous(
        realtype t, N_Vector u_and_s, N_Vector dydt_vec, void *user_data) {
  
  realtype* u              = N_VGetArrayPointer(u_and_s);
  realtype* s              = u + NEQ;
  realtype* dydt           = N_VGetArrayPointer(dydt_vec);
  realtype* dsdt           = dydt + NEQ;
  realtype* parameters     = ((UserData) user_data)->parameters;
  
  dydt_ode(u,dydt,parameters);
  
  double cx = DX * (N_MESH_POINTS + 1) * (N_MESH_POINTS+1) / (L*L);
  double cy = DY * (N_MESH_POINTS + 1) * (N_MESH_POINTS+1) / (L*L);
  
  double jac_xx = - 2 * cx - (B+1) + 2 * X(0) * Y(0);
  double jac_xy = X(0) * X(0);
  DSX_DT(0)     = jac_xx * SX(0) + cx * SX(1) + jac_xy * SY(0);
  
	double jac_yy = -2 * cy - X(0) * X(0);
  double jac_yx = B - 2 * X(0) * Y(0);
	DSY_DT(0)     = jac_yy * SY(0) + cy * SY(1) + jac_yx * SX(0);
  
  for (int i = 1; i < N_MESH_POINTS - 1; i++) {
    jac_xx    = - 2 * cx     - (B+1) + 2 * X(i) * Y(i);
    jac_xy    = X(i) * X(i);
    DSX_DT(i) = cx * SX(i-1) + jac_xx * SX(i) + cx * SX(i+1) + jac_xy * SY(i);
    
    jac_yy    = - 2 * cy         - X(i) * X(i);
    jac_yx    = B - 2 * X(i) * Y(i);
    DSY_DT(i) = cy * SY(i-1) + jac_yy * SY(i) + cy * SY(i+1) + jac_yx * SX(i);
  }
  
  int i = N_MESH_POINTS - 1;
  
  jac_xx    = - 2 * cx     - (B+1) + 2 * X(i) * Y(i);
  jac_xy    = X(i) * X(i);
  DSX_DT(i) = cx * SX(i-1) + jac_xx * SX(i) + jac_xy * SY(i);
  
  jac_yy    = - 2 * cy         - X(i) * X(i);
  jac_yx    = B - 2 * X(i) * Y(i);
  DSY_DT(i) = cy * SY(i-1) + jac_yy * SY(i) + jac_yx * SX(i);
  
  return 0;
}

#define X_INDEX(i) (2*(i))
#define Y_INDEX(i) (2*(i)+1)

#if DENSE_JACOBIAN
#define JAC_XX(i,j) SM_ELEMENT_D(jacobian, X_INDEX(i), X_INDEX(j))
#define JAC_XY(i,j) SM_ELEMENT_D(jacobian, X_INDEX(i), Y_INDEX(j))
#define JAC_YX(i,j) SM_ELEMENT_D(jacobian, Y_INDEX(i), X_INDEX(j))
#define JAC_YY(i,j) SM_ELEMENT_D(jacobian, Y_INDEX(i), Y_INDEX(j))
#endif

#if BANDED_JACOBIAN
#define JAC_XX(i,j) SM_ELEMENT_B(jacobian, X_INDEX(i), X_INDEX(j))
#define JAC_XY(i,j) SM_ELEMENT_B(jacobian, X_INDEX(i), Y_INDEX(j))
#define JAC_YX(i,j) SM_ELEMENT_B(jacobian, Y_INDEX(i), X_INDEX(j))
#define JAC_YY(i,j) SM_ELEMENT_B(jacobian, Y_INDEX(i), Y_INDEX(j))
#endif

#define NMP N_MESH_POINTS

static int jacobian_simultaneous(realtype t, N_Vector u_and_s, N_Vector fy,
              SUNMatrix jacobian, void* user_data,
              N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
  
  realtype* u   = N_VGetArrayPointer(u_and_s);
  realtype* parameters = ((UserData) user_data)->parameters;
  
  double cx = DX * (N_MESH_POINTS + 1) * (N_MESH_POINTS+1) / (L*L);
  double cy = DY * (N_MESH_POINTS + 1) * (N_MESH_POINTS+1) / (L*L);
   #if DENSE_JACOBIAN
  if (SUNMatGetID(jacobian) != SUNMATRIX_DENSE) {
    mexErrMsgIdAndTxt("cvode:wrong_matrix_type",
            "expected type %d = SUNMATRIX_DENSE, got type %d",
            SUNMATRIX_DENSE, SUNMatGetID(jacobian));
  }
  #endif
  
  #if BANDED_JACOBIAN
  if (SUNMatGetID(jacobian) != SUNMATRIX_BAND) {
    mexErrMsgIdAndTxt("cvode:wrong_matrix_type",
            "expected type %d = SUNMATRIX_BAND, got type %d",
            SUNMATRIX_DENSE, SUNMatGetID(jacobian));
  }
  #endif
  
  for ( int i = 0; i < N_MESH_POINTS; i++ ) {
    JAC_XX(i,i) = - 2 * cx - (B+1) + 2 * X(i) * Y(i);
    JAC_XY(i,i) =   X(i) * X(i);
    JAC_YX(i,i) =   B - 2 * X(i) * Y(i);
    JAC_YY(i,i) = - 2 * cy - X(i) * X(i);
  }
  
  for ( int i = 0; i < N_MESH_POINTS - 1; i++ ) {
    JAC_XX(i    , i + 1) = cx;
    JAC_XX(i + 1, i    ) = cx;
    JAC_YY(i    , i + 1) = cy;
    JAC_YY(i + 1, i    ) = cy;
  }
   
  for ( int i = 0; i < N_MESH_POINTS; i++ ) {
    JAC_XX(NMP + i, NMP + i) = - 2 * cx - (B+1) + 2 * X(i) * Y(i);
    JAC_XY(NMP + i, NMP + i) =   X(i) * X(i);
    JAC_YX(NMP + i, NMP + i) =    B - 2 * X(i) * Y(i);
    JAC_YY(NMP + i, NMP + i) = - 2 * cy - X(i) * X(i);
  }
  
  for ( int i = 0; i < N_MESH_POINTS - 1; i++ ) {
    JAC_XX(NMP + i    , NMP + i + 1) = cx;
    JAC_XX(NMP + i + 1, NMP + i    ) = cx;
    JAC_YY(NMP + i    , NMP + i + 1) = cy;
    JAC_YY(NMP + i + 1, NMP + i    ) = cy;
  }
  
  return 0;
}
static void error_handler(int error_code, const char* module,
         const char* function, char* message, void *eh_data) {
 mexPrintf("Error in %s() in sundails module %s. error code: %d: %s.\n  %s\n",
           function, module,
           error_code, CVodeGetReturnFlagName(error_code), message);
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

static void PrintFinalStats(void *cvode, booleantype sensi,
                            booleantype err_con, int sensi_meth)
{
  long int nst;
  long int nfe, nsetups, nni, ncfn, netf;
  long int nfSe, nfeS, nsetupsS, nniS, ncfnS, netfS;
  int flag;

  flag = CVodeGetNumSteps(cvode, &nst);
  //check_flag(&flag, "CVodeGetNumSteps", 1);
  flag = CVodeGetNumRhsEvals(cvode, &nfe);
  //check_flag(&flag, "CVodeGetNumRhsEvals", 1);
  flag = CVodeGetNumLinSolvSetups(cvode, &nsetups);
  //check_flag(&flag, "CVodeGetNumLinSolvSetups", 1);
  flag = CVodeGetNumErrTestFails(cvode, &netf);
  //check_flag(&flag, "CVodeGetNumErrTestFails", 1);
  flag = CVodeGetNumNonlinSolvIters(cvode, &nni);
  //check_flag(&flag, "CVodeGetNumNonlinSolvIters", 1);
  flag = CVodeGetNumNonlinSolvConvFails(cvode, &ncfn);
  //check_flag(&flag, "CVodeGetNumNonlinSolvConvFails", 1);

  if (sensi) {
    flag = CVodeGetSensNumRhsEvals(cvode, &nfSe);
    //check_flag(&flag, "CVodeGetSensNumRhsEvals", 1);
    flag = CVodeGetNumRhsEvalsSens(cvode, &nfeS);
    //check_flag(&flag, "CVodeGetNumRhsEvalsSens", 1);
    flag = CVodeGetSensNumLinSolvSetups(cvode, &nsetupsS);
    //check_flag(&flag, "CVodeGetSensNumLinSolvSetups", 1);
    if (err_con) {
      flag = CVodeGetSensNumErrTestFails(cvode, &netfS);
      //check_flag(&flag, "CVodeGetSensNumErrTestFails", 1);
    } else {
      netfS = 0;
    }
    if ((sensi_meth == CV_STAGGERED) || (sensi_meth == CV_STAGGERED1)) {
      flag = CVodeGetSensNumNonlinSolvIters(cvode, &nniS);
      //check_flag(&flag, "CVodeGetSensNumNonlinSolvIters", 1);
      flag = CVodeGetSensNumNonlinSolvConvFails(cvode, &ncfnS);
      //check_flag(&flag, "CVodeGetSensNumNonlinSolvConvFails", 1);
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



