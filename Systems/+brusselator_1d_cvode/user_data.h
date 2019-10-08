#define N_PARAMETERS          5
#define N_MESH_POINTS         30
#define NEQ                   2 * N_MESH_POINTS

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