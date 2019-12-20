function recompile_for_octave(N_MESH_POINTS)
  cl_matcontL_path = System_of_ODEs.get_cl_matcontL_path();
  cvodes_path = fullfile(cl_matcontL_path, 'SystemFileGenerator', 'cvodes');
  system_path = get_path();
  generator_path = fullfile(cl_matcontL_path, 'SystemFileGenerator')
  
  

  mex( ... 
    sprintf('-DN_MESH_POINTS=%d', N_MESH_POINTS), ...
    fullfile(generator_path, 'templates', 'cvode_mex.c'), ... 
    fullfile(system_path, 'dydt_ode.c'), ... 
    fullfile(system_path, 'dydt_cvode.c'), ... 
    fullfile(system_path, 'jacobian_cvode.c'), ... 
    fullfile(system_path, 'd_sensitivity_dt.c'), ... 
    fullfile(system_path, 'd_sensitivity_dt_pars.c'), ... 
    ['-I' fullfile(cvodes_path, 'include/')], ... 
    ['-I' system_path '/'],  ... 
    '-v', ...
    '-Wa,-g3', ... 
    [ cvodes_path '/cvodes/cvodea.c'], ... 
    [ cvodes_path '/cvodes/cvodea_io.c'], ... 
    [ cvodes_path '/cvodes/cvodes.c'], ... 
    [ cvodes_path '/cvodes/cvodes_bandpre.c'], ... 
    [ cvodes_path '/cvodes/cvodes_bbdpre.c'], ... 
    [ cvodes_path '/cvodes/cvodes_diag.c'], ... 
    [ cvodes_path '/cvodes/cvodes_direct.c'], ... 
    [ cvodes_path '/cvodes/cvodes_io.c'], ... 
    [ cvodes_path '/cvodes/cvodes_ls.c'], ... 
    [ cvodes_path '/cvodes/cvodes_nls.c'], ... 
    [ cvodes_path '/cvodes/cvodes_nls_sim.c'], ... 
    [ cvodes_path '/cvodes/cvodes_nls_stg.c'], ... 
    [ cvodes_path '/cvodes/cvodes_nls_stg1.c'], ... 
    [ cvodes_path '/cvodes/cvodes_spils.c'], ... 
    [ cvodes_path '/nvector/serial/fnvector_serial.c'], ... 
    [ cvodes_path '/nvector/serial/nvector_serial.c'], ... 
    [ cvodes_path '/sundials/sundials_band.c'], ... 
    [ cvodes_path '/sundials/sundials_dense.c'], ... 
    [ cvodes_path '/sundials/sundials_direct.c'], ... 
    [ cvodes_path '/sundials/sundials_iterative.c'], ... 
    [ cvodes_path '/sundials/sundials_linearsolver.c'], ... 
    [ cvodes_path '/sundials/sundials_math.c'], ... 
    [ cvodes_path '/sundials/sundials_matrix.c'], ... 
    [ cvodes_path '/sundials/sundials_nonlinearsolver.c'], ... 
    [ cvodes_path '/sundials/sundials_nvector.c'], ... 
    [ cvodes_path '/sundials/sundials_nvector_senswrapper.c'], ... 
    [ cvodes_path '/sundials/sundials_version.c'], ... 
    [ cvodes_path '/sunlinsol/band/fsunlinsol_band.c'], ... 
    [ cvodes_path '/sunlinsol/band/sunlinsol_band.c'], ... 
    [ cvodes_path '/sunlinsol/dense/fsunlinsol_dense.c'], ... 
    [ cvodes_path '/sunlinsol/dense/sunlinsol_dense.c'], ... 
    [ cvodes_path '/sunmatrix/band/fsunmatrix_band.c'], ... 
    [ cvodes_path '/sunmatrix/band/sunmatrix_band.c'], ... 
    [ cvodes_path '/sunmatrix/dense/fsunmatrix_dense.c'], ... 
    [ cvodes_path '/sunmatrix/dense/sunmatrix_dense.c'], ... 
    [ cvodes_path '/sunmatrix/sparse/fsunmatrix_sparse.c'], ... 
    [ cvodes_path '/sunmatrix/sparse/sunmatrix_sparse.c'], ... 
    [ cvodes_path '/sunnonlinsol/fixedpoint/fsunnonlinsol_fixedpoint.c'], ... 
    [ cvodes_path '/sunnonlinsol/fixedpoint/sunnonlinsol_fixedpoint.c'], ... 
    [ cvodes_path '/sunnonlinsol/newton/fsunnonlinsol_newton.c'], ... 
    [ cvodes_path '/sunnonlinsol/newton/sunnonlinsol_newton.c'], ... 
    '-output', ... 
    fullfile(system_path, 'cvode') ... 
 )

mex( ...
    fullfile(system_path, 'dydt_mex.c'), ...
    fullfile(system_path, 'dydt_ode.c'), ...
    sprintf('-DN_MESH_POINTS=%d', N_MESH_POINTS), ...
    '-output', ...
    fullfile(system_path, 'dydt_mex'));
  
mex( ...
    fullfile(system_path, 'jacobian_mex.c'), ...
    sprintf('-DN_MESH_POINTS=%d', N_MESH_POINTS), ...
    '-output', ...
    fullfile(system_path, 'jacobian_mex'));
  