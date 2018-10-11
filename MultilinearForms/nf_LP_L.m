function a = nf_LP_L(CISdata,x,p,borders)  % MP
%
% nf_LP_L(A0, x, p)
% compute fold normal form coefficient at a fold bifurcation
% A0 q0 = 0      A0^T p0 = 0
% coef = < p0, B(q0, q0) >
%

global contopts 

A0   = CISdata.A;
T    = CISdata.T;
NSub = CISdata.NSub;
Asub   = T(1:NSub,1:NSub); 
Q0     = CISdata.Q(:,1:NSub);

if nargin < 4 % compute borders
    [V, D] = eig(Asub);
    [~, ind] = min(abs(diag(D)));
    borders.v = V(:, ind);
    
    [W, D] = eig(Asub');
    [~, ind] = min(abs(diag(D)));
    borders.w = W(:, ind);
end

Bord    = [Asub borders.w;borders.v' 0];
bunit   = [zeros(NSub,1); 1];
vext    = Bord\bunit;
v1      = vext(1:end-1);
q0      = Q0 * v1;

Bord2    = [A0' Q0*borders.v; borders.w'* Q0' 0];
bunit2   = [zeros(length(A0),1); 1];
Wext     = Bord2\bunit2;
p0       = Wext(1:end-1);

% normalize 
q0=q0/norm(q0);                 % <q0, q0> = 1
p0=p0/(q0'*p0);                 % <q0, p0> = 1

% normal form coefficient
hessIncrement = (contopts.Increment)^(3.0/4.0);
Bq0q0 = multilinear2(q0,q0,x,p,hessIncrement);	% B(q0,q0)
a = p0'*Bq0q0/2.0;