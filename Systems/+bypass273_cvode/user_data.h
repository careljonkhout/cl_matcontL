#define N_PARAMETERS          1
#define NEQ                   273
#define ANALYTIC_JACOBIAN     true

#define DENSE   0
#define BANDED  1
#define SPARSE  2

#define JACOBIAN_STORAGE      DENSE

#define UPPER_BANDWIDTH       265
#define LOWER_BANDWIDTH       267

typedef struct {
  double* parameters;
} *UserData;
