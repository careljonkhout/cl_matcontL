

CC       = /usr/bin/gcc
DEBUG    = -O3
INCLUDES = -I"/home/carel/Documents/sundails/instdir/include" -I.  -I"/usr/local/MATLAB/R2019a/extern/include" -I"/usr/local/MATLAB/R2019a/simulink/include"
CFLAGS   = -c -DMATLAB_DEFAULT_RELEASE=R2017b -DUSE_MEX_CMD  -D_GNU_SOURCE -DMATLAB_MEX_FILE  -fexceptions -fPIC -fno-omit-frame-pointer -pthread $(DEBUG) $(INCLUDES)

default:  cvode_mex.o dydt_cvode.o jacobian_cvode.o d_sensitivity_dt.o c_mexapi_version.o
	$(CC) -pthread -Wl,--no-undefined -Wl,-rpath-link,/usr/local/MATLAB/R2019a/bin/glnxa64 -shared  -g -Wl,--version-script,"/usr/local/MATLAB/R2019a/extern/lib/glnxa64/c_exportsmexfileversion.map" object_files/cvode_mex.o object_files/dydt_cvode.o object_files/jacobian_cvode.o object_files/d_sensitivity_dt.o object_files/c_mexapi_version.o   /home/carel/Documents/sundails/builddir/src/cvodes/libsundials_cvodes.a /home/carel/Documents/sundails/builddir/src/nvector/serial/libsundials_nvecserial.a /home/carel/Documents/sundails/builddir/src/nvector/openmp/libsundials_nvecopenmp.a /home/carel/Documents/sundails/builddir/src/sunlinsol/band/libsundials_sunlinsolband.a /home/carel/Documents/sundails/builddir/src/sunlinsol/dense/libsundials_sunlinsoldense.a /home/carel/Documents/sundails/builddir/src/sunlinsol/pcg/libsundials_sunlinsolpcg.a /home/carel/Documents/sundails/builddir/src/sunlinsol/spbcgs/libsundials_sunlinsolspbcgs.a /home/carel/Documents/sundails/builddir/src/sunlinsol/spfgmr/libsundials_sunlinsolspfgmr.a /home/carel/Documents/sundails/builddir/src/sunlinsol/spgmr/libsundials_sunlinsolspgmr.a /home/carel/Documents/sundails/builddir/src/sunlinsol/sptfqmr/libsundials_sunlinsolsptfqmr.a /home/carel/Documents/sundails/builddir/src/sunnonlinsol/fixedpoint/libsundials_sunnonlinsolfixedpoint.a /home/carel/Documents/sundails/builddir/src/sunnonlinsol/newton/libsundials_sunnonlinsolnewton.a /home/carel/Documents/sundails/builddir/src/sunmatrix/band/libsundials_sunmatrixband.a /home/carel/Documents/sundails/builddir/src/sunmatrix/dense/libsundials_sunmatrixdense.a /home/carel/Documents/sundails/builddir/src/sunmatrix/sparse/libsundials_sunmatrixsparse.a  -L"/usr/local/MATLAB/R2019a/bin/glnxa64" -fopenmp -lmx -lmex -lmat -lm -lstdc++ -o /home/carel/Documents/cl_matcontL/SystemFileGenerator/../Systems/+brusselator_1d_N_10/cvode.mexa64

cvode_mex.o: object_files
	$(CC) $(CFLAGS) ../../SystemFileGenerator/cvode/cvode_mex.c -o ./object_files/cvode_mex.o
dydt_cvode.o: object_files
	$(CC) $(CFLAGS) dydt_cvode.c       -o object_files/dydt_cvode.o
jacobian_cvode.o: object_files
	$(CC) $(CFLAGS) jacobian_cvode.c   -o object_files/jacobian_cvode.o
d_sensitivity_dt.o: object_files
	$(CC) $(CFLAGS) d_sensitivity_dt.c -o object_files/d_sensitivity_dt.o	
c_mexapi_version.o: object_files
	$(CC) $(CFLAGS) "/usr/local/MATLAB/R2019a/extern/version/c_mexapi_version.c" -o object_files/c_mexapi_version.o

object_files:
	mkdir object_files
	

