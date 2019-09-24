#define N_PARAMETERS          5
#define N_MESH_POINTS         50
#define NEQ                   2*N_MESH_POINTS
#define ANALYTIC_JACOBIAN     true
#define BANDED_JACOBIAN       true
#define DENSE_JACOBIAN        !BANDED_JACOBIAN
#define CVODES                1
#define MANUAL                2
#define STAGGERED             3
#define DEFAULT_SENSITIVITY   CVODES

typedef struct {
  double* parameters;
  double temp[NEQ];
} *UserData;