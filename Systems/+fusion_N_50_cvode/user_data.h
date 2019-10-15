#define N_PARAMETERS 3
#define NEQ          147
#define ANALYTIC_JACOBIAN false
#define BANDED_JACOBIAN       true
#define DENSE_JACOBIAN        !BANDED_JACOBIAN
#define UPPER_HALF_BANDWIDTH  10
#define LOWER_HALF_BANDWIDTH  10


typedef struct {
  double* parameters;
} *UserData;
