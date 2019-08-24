#define N_PARAMETERS          5
#define N_MESH_POINTS         200
#define NEQ                   2*N_MESH_POINTS
#define ANALYTIC_JACOBIAN     true
#define BANDED_JACOBIAN       true
#define DENSE_JACOBIAN        !BANDED_JACOBIAN
#define MANUAL_SENSITIVITY    true

typedef struct {
  double* parameters;
} *UserData;