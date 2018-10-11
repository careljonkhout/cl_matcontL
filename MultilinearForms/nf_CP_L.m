function c = nf_CP_L(CISdata,x,p)
%
% compute cusp normal form coefficient.
%
global cds contopts

NSub   = CISdata.NSub;
A0     = CISdata.A;
T      = CISdata.T;
Asub   = T(1:NSub,1:NSub);
Q0     = CISdata.Q(:,1:NSub);

if isequal(cds.curve, @limitpointL)
    borders = cds.borders; % DV: use borders from LP curve
else
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

% compute coefficients
hessIncrement = (contopts.Increment)^(3.0/4.0);
ten3Increment = (contopts.Increment)^(3.0/5.0);

Bq0q0 = multilinear2(q0,q0,x,p,hessIncrement);
Cq0q0q0 = multilinear3(q0,q0,q0,x,p,ten3Increment);

Bord = [A0 p0 ; q0' 0];
a  = p0'*Bq0q0/2.0;
h2 = Bord\[-Bq0q0 + 2*a*q0; 0];%-A^{INV}(B(q0,q0)-<p0,B(q0,q0)>q0)
h2 = h2(1:end-1);
Bq0h2 = multilinear2(q0,h2,x,p,hessIncrement);

c = p0'*(Cq0q0q0 + 3*Bq0h2)/6.0;