#include <sundials/sundials_types.h>
#include <string.h>

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

void check_double(Argument arg) {
  if ( !mxIsDouble(arg.value) ) {
    mexErrMsgIdAndTxt("mex_utils:not_double", 
             "%s is not a vector of doubles.", arg.name);
  }
}



void check_size(Argument arg, int size) {
  if ( mxGetNumberOfElements(arg.value) != size ) {
    mexErrMsgIdAndTxt("mex_utils:incorrect_size",  
            "Input vector %s does not have the correct size. "
            "The correct size is %d. The actual size is %d.",
            arg.name, size, (int)mxGetNumberOfElements(arg.value));
  }
}

void check_scalar(Argument arg) {
  if ( mxGetNumberOfElements(arg.value) != 1 ) {
    mexErrMsgIdAndTxt("mex_utils:not_scalar",  
            "%s is not scalar, but it is an array of size %d.",
            arg.name, (int)mxGetNumberOfElements(arg.value));
  }
}

void check_real(Argument arg) {
  if ( mxIsComplex(arg.value) ) {
    mexErrMsgIdAndTxt("mex_utils:not_real", "%s is not real.", arg.name);
  }
}

void check_positive(Argument arg) {
  if (mxGetScalar(arg.value) <= 0) {
    mexErrMsgIdAndTxt("mex_utils:not_positive",  
            "Input %s is not positive.", arg.name);
  }
}

double* my_mex_get_doubles(const mxArray* array) {
  #if MX_HAS_INTERLEAVED_COMPLEX
  return (double*) mxGetDoubles(array);
  #else
  return (double*) mxGetPr(     array);
  #endif 
}

double* get_doubles(Argument arg, int size) {
  check_double(arg);
  check_real  (arg);
  check_size  (arg, size);
  return my_mex_get_doubles(arg.value);
}

double get_scalar(Argument arg) {
  check_double(arg);
  check_real  (arg);
  check_scalar(arg);
  return mxGetScalar(arg.value);
}

double get_positive_scalar(Argument arg) {
  check_double  (arg);
  check_real    (arg);
  check_scalar  (arg);
  check_positive(arg);
  return mxGetScalar(arg.value);
}

int get_nonnegative_int(Argument arg) {
  double value = get_scalar(arg);
  if (round(value) != value) {
    mexErrMsgIdAndTxt("mex_utils:not_integer",  
            "Input %s is not a whole number.", arg.name);
  }
  if (value < 0) {
    mexErrMsgIdAndTxt("mex_utils:negative",  
            "Input %s is negative.", arg.name);
  }
  // note: casting to int does not round, but truncates
  // therefore we add 0.5
  return (int) (value + 0.5);
}

booleantype get_bool(Argument arg) {
  if ( !mxIsLogicalScalar(arg.value) ) {
    mexErrMsgIdAndTxt("mex_utils:not_bool", 
      "Input %s is not a boolean.", arg.name);
  }
  return mxIsLogicalScalarTrue(arg.value);
}