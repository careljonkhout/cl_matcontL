#include <sundials/sundials_types.h>
#include <string.h>
#include <mex.h>
#include <math.h>

typedef struct {
  char* name;
  const mxArray* value;
} Argument;


void check_n_inputs(int n_inputs, int required, char* function) {
  if (n_inputs != required) {
    mexErrMsgIdAndTxt("mex_utils:n_inputs",
            "%s requires %d inputs, but got %d inputs",
            function, required, n_inputs);
  }
}

void check_double(const mxArray* array, char* arrayname) {
  if ( !mxIsDouble(array) ) {
    mexErrMsgIdAndTxt("mex_utils:not_double", 
             "%s is not a vector of doubles.", arrayname);
  }
}

void check_bool(const mxArray* array, char* arrayname) {
  if ( !mxIsLogicalScalar(array) ) {
    mexErrMsgIdAndTxt("mex_utils:not_bool", 
            "Input %s is not a boolean.", arrayname);
  }
}

void check_size(const mxArray* array, int size, char* arrayname) {
  if ( mxGetNumberOfElements(array) != size ) {
    mexErrMsgIdAndTxt("mex_utils:incorrect_size",  
            "Input vector %s does not have the correct size. "
            "The correct size is %d. The actual size is %d.",
            arrayname, size, (int)mxGetNumberOfElements(array));
  }
}

void check_scalar(const mxArray* array, char* arrayname) {
  if ( mxGetNumberOfElements(array) != 1 ) {
    mexErrMsgIdAndTxt("mex_utils:not_scalar",  
            "%s is not scalar, but it is an array of size %d.",
            arrayname, (int)mxGetNumberOfElements(array));
  }
}

void check_real(const mxArray* array, char* arrayname) {
  if ( mxIsComplex(array) ) {
    mexErrMsgIdAndTxt("mex_utils:not_real", "%s is not real.", arrayname);
  }
}

double* get_parameters(const mxArray* inputs[],
        double* parameters, int n_parameters) {
  for ( int i = 0; i < n_parameters; i++ ) {
    check_double(inputs[i + 2], "One of the parameters");
    check_real  (inputs[i + 2], "One of the parameters");
    check_scalar(inputs[i + 2], "One of the parameters");
    parameters[i] = mxGetScalar(inputs[i+2]);
  }
}

void check_positive(const mxArray* scalar, char* scalarname) {
  if (mxGetScalar(scalar) <= 0) {
    mexErrMsgIdAndTxt("mex_utils:not_positive",  
            "Input %s is not positive.", scalarname);
  }
}

double* my_mex_get_doubles(const mxArray* array) {
  #if MX_HAS_INTERLEAVED_COMPLEX
  return (double*) mxGetDoubles(array);
  #else
  return (double*) mxGetPr(     array);
  #endif 
}

double* get_doubles(const mxArray* array, char* input_name, int size) {
  check_double(array,       input_name);
  check_real  (array,       input_name);
  check_size  (array, size, input_name);
  return my_mex_get_doubles(array);
}

double get_scalar(const mxArray* array, char* input_name) {
  check_double(array, input_name);
  check_real  (array, input_name);
  check_scalar(array, input_name);
  return mxGetScalar(array);
}


double get_positive_scalar(const mxArray* array, char* input_name) {
  check_double  (array, input_name);
  check_real    (array, input_name);
  check_scalar  (array, input_name);
  check_positive(array, input_name);
  return mxGetScalar(array);
}

int get_int(const mxArray* array, char* input_name) {
  double value = get_scalar(array, input_name);
  if (round(value) != value) {
    mexErrMsgIdAndTxt("mex_utils:not_integer",  
            "Input %s is not a whole number.", input_name);
  }
  if (value < 0) {
    mexErrMsgIdAndTxt("mex_utils:negative",  
            "Input %s is negative.", input_name);
  }
  // note: casting to int does not round, but truncates
  // therefore we add 0.5
  return (int) (value + 0.5);
}


booleantype get_positive_scalar_if_name_matches(Argument arg, char* name,
                                                              realtype* ptr) {
  booleantype matches = ! strcmp(name, arg.name);
  if ( matches ) {
    *ptr = (realtype) get_positive_scalar(arg.value, arg.name);
  }
  return matches;
}



