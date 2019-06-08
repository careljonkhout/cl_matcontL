function [x0,v0]= init_H_H_L(probfile, u, p, ap, data)
%
% Initializes a Hopf bifurcation continuation from a Hopf point found 
%
global cds

% check input
if nargin > 5
    error('init_H_H_L cannot be called with more than 5 arguments');
end
if isempty(data)
    error('Data from equilibriumL must be present'); 
end

if isempty(u)
    u = data.x(1:end-1);
end
if isempty(p)
    p = data.P0;
    p(data.ap) = data.x(end); % DV 2018
    if size(p, 2) > 1
        p = p';
    end
end
if size(ap,2)~=2
    errordlg('Two active parameter are needed for a Hopfpoint bifurcation continuation');
end

x0 = [u; p(ap); data.s.kapa];        % change continuation parameters
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
cds.ncoo         = length(x0) - 3;
cds.nap          = length(ap);
cds.usernorm     = user_norm;
cds.userf        = user_func;
cds.curve        = @hopfL;