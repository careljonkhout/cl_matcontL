#define N_PARAMETERS          7
#define NEQ                   100
#define ANALYTIC_JACOBIAN     true

#define DENSE   0
#define BANDED  1
#define SPARSE  2

#define JACOBIAN_STORAGE      BANDED

#define UPPER_BANDWIDTH       2
#define LOWER_BANDWIDTH       2

#define OPEN_MP               false
#define N_THREADS             2

typedef struct {
  double* parameters;
} *UserData;
