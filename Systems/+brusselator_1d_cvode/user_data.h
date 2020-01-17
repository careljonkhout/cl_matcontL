#define N_PARAMETERS          6
#define NEQ                   2 * N_MESH_POINTS

#define DENSE   0
#define BANDED  1
#define SPARSE  2

#define JACOBIAN_STORAGE      BANDED

#define ANALYTIC_JACOBIAN     true

#define UPPER_BANDWIDTH       2
#define LOWER_BANDWIDTH       2

#define OPEN_MP               false
#define N_THREADS             2

#define PRINT_TIME            true

#define N_MESH_POINTS_AS_PARAMETER true
#define PDE_DIMENSION              2

typedef struct {
  double* parameters;
} *UserData;