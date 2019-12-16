#define N_PARAMETERS          5
#define NEQ                   5
#define ANALYTIC_JACOBIAN     true

#define DENSE   0
#define BANDED  1
#define SPARSE  2

#define JACOBIAN_STORAGE      DENSE

#define UPPER_BANDWIDTH       3
#define LOWER_BANDWIDTH       1

#define OPEN_MP               false
#define N_THREADS             2

typedef struct {
  realtype* parameters;
} *UserData;
