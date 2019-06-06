N=25;
odefile = str2func(sprintf('fusion_N_%d_hess', N));
handles = feval(odefile);
jacobian = handles{3};

x = ones(3*(N-1),1);
p = {-1,-0.3,-0.7};
ap = 3;
jac1 = feval(jacobian,0,x,p{:});

odefile2 = str2func(sprintf('fusion_precomputed_with_sage_N_%d', N));
handles2 = feval(odefile2);
jac2 = feval(handles2{3},0,x,p{:});
spy(abs(jac1-jac2) < 1e-13)
