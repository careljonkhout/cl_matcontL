#define N_PARAMETERS 1
#define NEQ          273
#define ANALYTIC_JACOBIAN false
#define BANDED_JACOBIAN       false
#define DENSE_JACOBIAN        !BANDED_JACOBIAN

typedef struct {
  double* parameters;
} *UserData;
