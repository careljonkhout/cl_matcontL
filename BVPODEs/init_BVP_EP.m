function [x0, v0] = init_BVP_EP(probfile, x0, p, ap)

handles = feval(probfile);

init_prob = handles{1};
func_prob = handles{2};
jacu_prob = handles{3};
jacp_prob = handles{4};
user_norm = handles{5};
for ii = 6:length(handles)
    user_func{ii-5} = handles{ii};
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

global cds
cds.probfile     = probfile;
cds.func         = func_prob;
cds.Jacobian     = jacu_prob;
cds.JacobianP    = jacp_prob;
cds.P0           = p;
cds.ActiveParams = ap;
cds.ndim         = length(x0);
cds.ncoo         = length(x0) - 1;
cds.nap          = length(ap);
cds.usernorm     = user_norm;
cds.userf        = user_func;