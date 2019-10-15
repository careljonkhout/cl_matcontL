#define N_PARAMETERS          3
#define NEQ                   72
#define ANALYTIC_JACOBIAN     true

#define DENSE   0
#define BANDED  1
#define SPARSE  2

#define JACOBIAN_STORAGE      BANDED

#define UPPER_BANDWIDTH       5
#define LOWER_BANDWIDTH       5

typedef struct {
  double* parameters;
} *UserData;
