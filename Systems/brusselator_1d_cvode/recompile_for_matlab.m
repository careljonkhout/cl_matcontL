function recompile_for_matlab(N_MESH_POINTS)
  cd(get_path())
  cd(fullfile('..','..'))

  mex( ... 
    sprintf('-DN_MESH_POINTS=%d', N_MESH_POINTS), ...
    'SystemFileGenerator/templates/cvode_mex.c', ... 
    'Systems/brusselator_1d_cvode/dydt_ode.c', ... 
    'Systems/brusselator_1d_cvode/dydt_cvode.c', ... 
    'Systems/brusselator_1d_cvode/jacobian_cvode.c', ... 
    'Systems/brusselator_1d_cvode/d_sensitivity_dt.c', ... 
    'Systems/brusselator_1d_cvode/d_sensitivity_dt_pars.c', ... 
    '-ISystemFileGenerator/cvodes/include', ... 
    '-ISystemFileGenerator/../Systems/brusselator_1d_cvode', ... 
    ...'-g', ...
    'SystemFileGenerator/cvodes/cvodes/cvodea.c', ... 
    'SystemFileGenerator/cvodes/cvodes/cvodea_io.c', ... 
    'SystemFileGenerator/cvodes/cvodes/cvodes.c', ... 
    'SystemFileGenerator/cvodes/cvodes/cvodes_bandpre.c', ... 
    'SystemFileGenerator/cvodes/cvodes/cvodes_bbdpre.c', ... 
    'SystemFileGenerator/cvodes/cvodes/cvodes_diag.c', ... 
    'SystemFileGenerator/cvodes/cvodes/cvodes_direct.c', ... 
    'SystemFileGenerator/cvodes/cvodes/cvodes_io.c', ... 
    'SystemFileGenerator/cvodes/cvodes/cvodes_ls.c', ... 
    'SystemFileGenerator/cvodes/cvodes/cvodes_nls.c', ... 
    'SystemFileGenerator/cvodes/cvodes/cvodes_nls_sim.c', ... 
    'SystemFileGenerator/cvodes/cvodes/cvodes_nls_stg.c', ... 
    'SystemFileGenerator/cvodes/cvodes/cvodes_nls_stg1.c', ... 
    'SystemFileGenerator/cvodes/cvodes/cvodes_spils.c', ... 
    'SystemFileGenerator/cvodes/nvector/serial/fnvector_serial.c', ... 
    'SystemFileGenerator/cvodes/nvector/serial/nvector_serial.c', ... 
    'SystemFileGenerator/cvodes/sundials/sundials_band.c', ... 
    'SystemFileGenerator/cvodes/sundials/sundials_dense.c', ... 
    'SystemFileGenerator/cvodes/sundials/sundials_direct.c', ... 
    'SystemFileGenerator/cvodes/sundials/sundials_iterative.c', ... 
    'SystemFileGenerator/cvodes/sundials/sundials_linearsolver.c', ... 
    'SystemFileGenerator/cvodes/sundials/sundials_math.c', ... 
    'SystemFileGenerator/cvodes/sundials/sundials_matrix.c', ... 
    'SystemFileGenerator/cvodes/sundials/sundials_nonlinearsolver.c', ... 
    'SystemFileGenerator/cvodes/sundials/sundials_nvector.c', ... 
    'SystemFileGenerator/cvodes/sundials/sundials_nvector_senswrapper.c', ... 
    'SystemFileGenerator/cvodes/sundials/sundials_version.c', ... 
    'SystemFileGenerator/cvodes/sunlinsol/band/fsunlinsol_band.c', ... 
    'SystemFileGenerator/cvodes/sunlinsol/band/sunlinsol_band.c', ... 
    'SystemFileGenerator/cvodes/sunlinsol/dense/fsunlinsol_dense.c', ... 
    'SystemFileGenerator/cvodes/sunlinsol/dense/sunlinsol_dense.c', ... 
    'SystemFileGenerator/cvodes/sunmatrix/band/fsunmatrix_band.c', ... 
    'SystemFileGenerator/cvodes/sunmatrix/band/sunmatrix_band.c', ... 
    'SystemFileGenerator/cvodes/sunmatrix/dense/fsunmatrix_dense.c', ... 
    'SystemFileGenerator/cvodes/sunmatrix/dense/sunmatrix_dense.c', ... 
    'SystemFileGenerator/cvodes/sunmatrix/sparse/fsunmatrix_sparse.c', ... 
    'SystemFileGenerator/cvodes/sunmatrix/sparse/sunmatrix_sparse.c', ... 
    'SystemFileGenerator/cvodes/sunnonlinsol/fixedpoint/fsunnonlinsol_fixedpoint.c', ... 
    'SystemFileGenerator/cvodes/sunnonlinsol/fixedpoint/sunnonlinsol_fixedpoint.c', ... 
    'SystemFileGenerator/cvodes/sunnonlinsol/newton/fsunnonlinsol_newton.c', ... 
    'SystemFileGenerator/cvodes/sunnonlinsol/newton/sunnonlinsol_newton.c', ... 
    '-output', ... 
    'Systems/brusselator_1d_cvode/cvode' ... 
  )

mex( ...
    'Systems/brusselator_1d_cvode/dydt_mex.c', ...
    'Systems/brusselator_1d_cvode/dydt_ode.c', ...
    sprintf('-DN_MESH_POINTS=%d', N_MESH_POINTS), ...
    '-output', ...
    'Systems/brusselator_1d_cvode/dydt_mex');
  
mex( ...
    'Systems/brusselator_1d_cvode/jacobian_mex.c', ...
    sprintf('-DN_MESH_POINTS=%d', N_MESH_POINTS), ...
    '-output', ...
    'Systems/brusselator_1d_cvode/jacobian_mex');
  