N=10;
odefile = str2func(sprintf('fusion_N_%d_hess', N));
handles = feval(odefile);
jacobian = handles{3};
global cds
cds.options.SymDerivativeP = 0;
x = ones(3*(N-1),1);
p = {-1,-0.3,-0.7};
ap = 3;
cds.nphase = 3*(N-1);
h_parameters_finite_differences = chessp(odefile,jacobian,[],x,p,ap);


h_finite_diffenreces = chess(odefile,jacobian,[],x,p,ap);

h_analytical = feval(handles{5},0, x, p{:});
