copyfile( ... 
  '/home/carel/Documents/cl_matcontL/SystemFileGenerator/cvode/cvode_mex.c',  ... 
  '/home/carel/Documents/cl_matcontL/SystemFileGenerator/../Systems/+brusselator_1d_N_500/cvode_mex.c' ... 
)
mex( ... 
  '/home/carel/Documents/cl_matcontL/SystemFileGenerator/../Systems/+brusselator_1d_N_500/cvode_mex.c', ... 
  '/home/carel/Documents/cl_matcontL/SystemFileGenerator/../Systems/+brusselator_1d_N_500/dydt_cvode.c', ... 
  '/home/carel/Documents/cl_matcontL/SystemFileGenerator/../Systems/+brusselator_1d_N_500/jacobian_cvode.c', ... 
  '/home/carel/Documents/cl_matcontL/SystemFileGenerator/../Systems/+brusselator_1d_N_500/d_sensitivity_dt.c', ... 
  '-g', ... 
  '/home/carel/Documents/sundails/builddir/src/cvodes/libsundials_cvodes.a', ... 
  '/home/carel/Documents/sundails/builddir/src/nvector/serial/libsundials_nvecserial.a', ... 
  '/home/carel/Documents/sundails/builddir/src/sunlinsol/band/libsundials_sunlinsolband.a', ... 
  '/home/carel/Documents/sundails/builddir/src/sunlinsol/dense/libsundials_sunlinsoldense.a', ... 
  '/home/carel/Documents/sundails/builddir/src/sunlinsol/pcg/libsundials_sunlinsolpcg.a', ... 
  '/home/carel/Documents/sundails/builddir/src/sunlinsol/spbcgs/libsundials_sunlinsolspbcgs.a', ... 
  '/home/carel/Documents/sundails/builddir/src/sunlinsol/spfgmr/libsundials_sunlinsolspfgmr.a', ... 
  '/home/carel/Documents/sundails/builddir/src/sunlinsol/spgmr/libsundials_sunlinsolspgmr.a', ... 
  '/home/carel/Documents/sundails/builddir/src/sunlinsol/sptfqmr/libsundials_sunlinsolsptfqmr.a', ... 
  '/home/carel/Documents/sundails/builddir/src/sunnonlinsol/fixedpoint/libsundials_sunnonlinsolfixedpoint.a', ... 
  '/home/carel/Documents/sundails/builddir/src/sunnonlinsol/newton/libsundials_sunnonlinsolnewton.a', ... 
  '/home/carel/Documents/sundails/builddir/src/sunmatrix/band/libsundials_sunmatrixband.a', ... 
  '/home/carel/Documents/sundails/builddir/src/sunmatrix/dense/libsundials_sunmatrixdense.a', ... 
  '/home/carel/Documents/sundails/builddir/src/sunmatrix/sparse/libsundials_sunmatrixsparse.a', ... 
  '-output', ... 
  '/home/carel/Documents/cl_matcontL/SystemFileGenerator/../Systems/+brusselator_1d_N_500/cvode' ... 
)