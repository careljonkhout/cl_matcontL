#define N_PARAMETERS          5
#define N_MESH_POINTS         50
#define NEQ                   2*N_MESH_POINTS
#define ANALYTIC_JACOBIAN     true
#define BANDED_JACOBIAN       true
#define DENSE_JACOBIAN        !BANDED_JACOBIAN
#define CVODES       0
#define MANUAL       1
#define STAGGERED    2
#define SENSITIVITY  MANUAL

typedef struct {
  double* parameters;
  double temp[NEQ];
} *UserData;