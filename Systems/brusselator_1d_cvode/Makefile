DEPS         = user_data.h
DEBUG        = -g
COMPILE_OPTS = -c -DMATLAB_DEFAULT_RELEASE=R2017b -DUSE_MEX_CMD  -D_GNU_SOURCE -DMATLAB_MEX_FILE -fexceptions -fPIC -fno-omit-frame-pointer -pthread -O3 $(DEBUG)
OPENMP       = -fopenmp

default: cvode_mex.o dydt_ode.o dydt_cvode.o jacobian_cvode.o d_sensitivity_dt.o c_mexapi_version.o
	/usr/bin/gcc -v -pthread -Wl,--no-undefined -Wl,-rpath-link,/usr/local/MATLAB/R2019a/bin/glnxa64 -shared  $(DEBUG) -Wl,--version-script,"/usr/local/MATLAB/R2019a/extern/lib/glnxa64/c_exportsmexfileversion.map" object_files/c_mexapi_version.o object_files/cvode_mex.o object_files/dydt_ode.o object_files/dydt_cvode.o object_files/jacobian_cvode.o object_files/d_sensitivity_dt.o    /home/carel/Documents/sundails/builddir/src/cvodes/libsundials_cvodes.a /home/carel/Documents/sundails/builddir/src/nvector/serial/libsundials_nvecserial.a /home/carel/Documents/sundails/builddir/src/nvector/openmp/libsundials_nvecopenmp.a /home/carel/Documents/sundails/builddir/src/sunlinsol/band/libsundials_sunlinsolband.a /home/carel/Documents/sundails/builddir/src/sunlinsol/dense/libsundials_sunlinsoldense.a /home/carel/Documents/sundails/builddir/src/sunlinsol/pcg/libsundials_sunlinsolpcg.a /home/carel/Documents/sundails/builddir/src/sunlinsol/spbcgs/libsundials_sunlinsolspbcgs.a /home/carel/Documents/sundails/builddir/src/sunlinsol/spfgmr/libsundials_sunlinsolspfgmr.a /home/carel/Documents/sundails/builddir/src/sunlinsol/spgmr/libsundials_sunlinsolspgmr.a /home/carel/Documents/sundails/builddir/src/sunlinsol/sptfqmr/libsundials_sunlinsolsptfqmr.a /home/carel/Documents/sundails/builddir/src/sunnonlinsol/fixedpoint/libsundials_sunnonlinsolfixedpoint.a /home/carel/Documents/sundails/builddir/src/sunnonlinsol/newton/libsundials_sunnonlinsolnewton.a /home/carel/Documents/sundails/builddir/src/sunmatrix/band/libsundials_sunmatrixband.a /home/carel/Documents/sundails/builddir/src/sunmatrix/dense/libsundials_sunmatrixdense.a /home/carel/Documents/sundails/builddir/src/sunmatrix/sparse/libsundials_sunmatrixsparse.a  -L"/usr/local/MATLAB/R2019a/bin/glnxa64" $(OPENMP) -lmx -lmex -lmat -lm -lstdc++ -o cvode.mexa64

cvode_mex.o:
	/usr/bin/gcc $(COMPILE_OPTS) $(OPENMP) -I"/home/carel/Documents/sundails/instdir/include" -I.  -I"/usr/local/MATLAB/R2019a/extern/include" -I"/usr/local/MATLAB/R2019a/simulink/include"  "/home/carel/Documents/cl_matcontL/SystemFileGenerator/cvode/cvode_mex.c" -o object_files/cvode_mex.o
dydt_cvode.o:
	/usr/bin/gcc $(COMPILE_OPTS) $(OPENMP) -I"/home/carel/Documents/sundails/instdir/include" -I.  -I"/usr/local/MATLAB/R2019a/extern/include" -I"/usr/local/MATLAB/R2019a/simulink/include" dydt_cvode.c -o object_files/dydt_cvode.o
dydt_ode.o:
	/usr/bin/gcc $(COMPILE_OPTS) $(OPENMP) -I"/home/carel/Documents/sundails/instdir/include" -I.  -I"/usr/local/MATLAB/R2019a/extern/include" -I"/usr/local/MATLAB/R2019a/simulink/include" dydt_ode.c -o object_files/dydt_ode.o
jacobian_cvode.o:
	/usr/bin/gcc $(COMPILE_OPTS) $(OPENMP) -I"/home/carel/Documents/sundails/instdir/include" -I.  -I"/usr/local/MATLAB/R2019a/extern/include" -I"/usr/local/MATLAB/R2019a/simulink/include" jacobian_cvode.c -o object_files/jacobian_cvode.o
d_sensitivity_dt.o:
	/usr/bin/gcc $(COMPILE_OPTS) $(OPENMP) -I"/home/carel/Documents/sundails/instdir/include" -I.  -I"/usr/local/MATLAB/R2019a/extern/include" -I"/usr/local/MATLAB/R2019a/simulink/include" d_sensitivity_dt.c -o object_files/d_sensitivity_dt.o
	
c_mexapi_version.o:
	/usr/bin/gcc $(COMPILE_OPTS) $(OPENMP) -I"/home/carel/Documents/sundails/instdir/include" -I.  -I"/usr/local/MATLAB/R2019a/extern/include" -I"/usr/local/MATLAB/R2019a/simulink/include" -fexceptions -fPIC -fno-omit-frame-pointer -pthread $(DEBUG) "/usr/local/MATLAB/R2019a/extern/version/c_mexapi_version.c" -o object_files/c_mexapi_version.o
	

