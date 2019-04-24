N=75;

odefile = str2func(sprintf('fusion_precomputed_with_sage_N_%d', N));
odefile_mex = str2func(sprintf('fusion_mex_N_%d', N));
odefile_mex_sparse = str2func(sprintf('fusion_mex_sparse_N_%d', N));

handles = feval(odefile);
handles_mex = feval(odefile_mex);
handles_mex_sparse = feval(odefile_mex_sparse);

dydt     = feval(handles    {2},0,ones(3*(N-1),1),1,1,1);
dydt_mex = feval(handles_mex{2},0,ones(3*(N-1),1),1,1,1);

disp(max(abs(dydt-dydt_mex)))


jac            = feval(handles           {3},0,ones(3*(N-1),1),1,1,1);
jac_mex        = feval(handles_mex       {3},0,ones(3*(N-1),1),1,1,1);
jac_mex_sparse = feval(handles_mex_sparse{3},0,ones(3*(N-1),1),1,1,1);

disp(max(max(abs(jac-jac_mex))))
disp(max(max(abs(jac-jac_mex_sparse))))
