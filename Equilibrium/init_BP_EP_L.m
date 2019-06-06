function [x0,v0,options]= init_BP_EP_L(probfile, x, p, ap, data)
%

% [x0,v0]= init_BP_EPL(odefile, x, p, s, h)
%
% Branch switching at an Branching Point (BP) detected on an equilibrium
% curve
%
% data from located branching point must be available
% If x, p or ap are empty the values from the data struct are used
%

if isempty(data)
    error('Data structure must be given')
end

if isempty(x)
    x = data.x(1:end-1);
end

if isempty(p)
    p = data.P0;
end

if isempty(ap)
    ap = data.ap;
end
if length(ap) ~= 1
    error('Requires only one active parameter')
end
if ap ~= data.ap
    error('Active parameter must be the same as on the equilibrium curve')
end
x0 = [x; p(ap)];

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

global cds
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
cds.curve        = @equilibriumL;


if isfield(data, 's') && isfield(data.s,'xnext')
    x0 = data.s.xnext;
    v0 = data.s.vnext;
elseif isfield(data.s, 'v2') || ~isempty(data.s.v2)
    v0 = data.s.v2;
else
    error('No tangent vector found in data. BP might be degenerate');
end
