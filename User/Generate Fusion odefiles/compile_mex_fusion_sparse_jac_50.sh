/usr/bin/gcc -c -DMATLAB_DEFAULT_RELEASE=R2017b  -DUSE_MEX_CMD   -D_GNU_SOURCE -DMATLAB_MEX_FILE  -I"/usr/local/MATLAB/R2018b/extern/include" -I"/usr/local/MATLAB/R2018b/simulink/include" -fexceptions -fPIC -fno-omit-frame-pointer -pthread -O3 -fwrapv -DNDEBUG "/home/carel/Documents/cl_matcontL/User/Generate Fusion odefiles/fusion_sparse_jacobian_N_50.c" -o fusion_sparse_jacobian_N_50.o
/usr/bin/gcc -c -DMATLAB_DEFAULT_RELEASE=R2017b  -DUSE_MEX_CMD   -D_GNU_SOURCE -DMATLAB_MEX_FILE  -I"/usr/local/MATLAB/R2018b/extern/include" -I"/usr/local/MATLAB/R2018b/simulink/include" -fexceptions -fPIC -fno-omit-frame-pointer -pthread -O3 -fwrapv -DNDEBUG "/usr/local/MATLAB/R2018b/extern/version/c_mexapi_version.c" -o c_mexapi_version.o
/usr/bin/gcc -pthread -Wl,--no-undefined -Wl,-rpath-link,/usr/local/MATLAB/R2018b/bin/glnxa64 -shared  -O -Wl,--version-script,"/usr/local/MATLAB/R2018b/extern/lib/glnxa64/c_exportsmexfileversion.map" fusion_sparse_jacobian_N_50.o c_mexapi_version.o   -L"/usr/local/MATLAB/R2018b/bin/glnxa64" -lmx -lmex -lmat -lm -lstdc++ -o fusion_sparse_jacobian_N_50.mexa64

