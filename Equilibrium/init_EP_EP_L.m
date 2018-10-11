function [x0, v0] = init_EP_EP_L(probfile, x0, p, ap)
%print_diag(5,'In init_EP_EP_L\n');

global cds
handles = feval(probfile);

init_prob = handles{1};
func_prob = handles{2};
jacu_prob = handles{3};
jacp_prob = handles{4};
hesu_prob = handles{5};
hesp_prob = handles{6};
der3_prob = handles{7};
der4_prob = handles{8};
der5_prob = handles{9};
user_norm = handles{10};
user_func = [];
for ii = 11:length(handles)
    user_func{ii-10} = handles{ii};
end

if isempty(x0)
    if ~isempty(init_prob)
        p2 = num2cell(p);
        x0 = feval(init_prob, p2{:});
    else
        error('No initial point found');
    end
end
x0 = [x0; p(ap)];
v0 = [];

% initialize cds
cds.probfile     = probfile;
cds.func         = func_prob;
cds.Jacobian     = jacu_prob;
cds.JacobianP    = jacp_prob;
cds.Hessian      = hesu_prob;
cds.HessianP     = hesp_prob;
cds.der3         = der3_prob;
cds.der4         = der4_prob;
cds.der5         = der5_prob;
cds.P0           = p;
cds.ActiveParams = ap;
cds.ndim         = length(x0);
cds.ncoo         = length(x0) - 1;
cds.nap          = length(ap);
cds.usernorm     = user_norm;
cds.userf        = user_func;