function [x0,v0]= init_LP_LP_L(probfile, u, p, ap, data)
%  Initializes a Limit Point / Fold continuation from a LP point

clear global
global cds

% check input
if nargin > 5
    error('init_LP_LP_L cannot be called with more than 5 arguments');
end
if isempty(u)
    if nargin < 5 % no data struct available
        error('No initial point found');
    else
        u = data.x(1:end-length(data.ap));
    end
end
if isempty(p)
    if nargin < 5 % no data struct available
        error('No parameters found');
    else
        p = data.P0;
    end
end
if size(ap,2)~=2
    errordlg('Two active parameter are needed for a Hopfpoint bifurcation continuation');
end
x0 = [u; p(ap)'];        % change continuation parameters
v0 = [];

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
cds.ncoo         = length(x0) - 2;
cds.nap          = length(ap);
cds.usernorm     = user_norm;
cds.userf        = user_func;
cds.curve        = @limitpointL;

% cds.lastwh       = data.lastwh;  DV: We do not need this anymore, because
% we can use the CIS algorithm to assure smooth borders
cds.BranchParams = [];   % DV: branch switching on LP curve not yet tested