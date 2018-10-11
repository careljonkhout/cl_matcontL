function [coef, kapa] = nf_H_L(CISdata,x,p,kapa,borders)
%
% nf_H_L(A0,x,p,k,vext,wext)
% calculates the first lyapunov coefficient
% 
global cds contopts

A0 = CISdata.A;
T  = CISdata.T;
NSub = CISdata.NSub;

T0     = T(1:NSub,1:NSub);
Q0     = CISdata.Q(:,1:NSub);

if nargin < 5 % compute borders and kapa
    [V, D] = eig(T0);
    d = diag(D);
    I     = find(imag(d) ~= 0);
    if isempty(I)
        print_diag(1, 'nf_H_L: All eigenvalues are real \n');
        coef = [];
        return
    end
    Rmin  = min(abs(real(d(I))));
    I1    = find(abs(real(d)) == Rmin);
    borders.v =  [real(V(:,I1(1))) imag(V(:,I1(1)))];
    
    [W, D] = eig(T0');
    d = diag(D);
    I     = find(imag(d) ~= 0);
    Rmin  = min(abs(real(d(I))));
    I2    = find(abs(real(d)) == Rmin);
    borders.w =  [real(W(:,I2(1))) imag(W(:,I2(1)))];
    
    kapa = abs(d(I2(1)))^2;
end

% right eigenvector
RED = T0*T0+kapa*eye(NSub);
Bord  = [RED borders.w; borders.v' zeros(2)];
bunit = [zeros(NSub,2);eye(2)];

vext = Bord\bunit;
v1   = vext(1:NSub,1);
v2   = vext(1:NSub,2);

omega=sqrt(kapa);
alpha=v1'*T0*v2 - 1i*omega*(v1'*v2);
beta=-v1'*T0*v1 + 1i*omega*(v1'*v1);
q0 = Q0 * (alpha*v1 + beta*v2);

% left eigenvector
RED2 = A0*A0+kapa*speye(length(A0));
Bord2 = [RED2' Q0*borders.v; borders.w'*Q0' zeros(2)];
bunit2 = [zeros(length(A0), 2); eye(2)];

wext = bordCIS1(Bord2, bunit2, 2);  
w1   = wext(1:cds.ncoo,1);
w2   = wext(1:cds.ncoo,2);

alpha=w1'*A0'*w2+1i*omega*(w1'*w2);
beta=-w1'*A0'*w1-1i*omega*(w1'*w1);
p0=alpha*w1+beta*w2;

q0=q0/norm(q0);
p0=p0/(q0'*p0);

% normal form coefficient
hessIncrement = (contopts.Increment)^(3.0/4.0);
ten3Increment = (contopts.Increment)^(3.0/5.0);

Bq0q0  = multilinear2(q0, q0 ,x,p,hessIncrement); 
Bq0cq0 = multilinear2(q0,conj(q0),x,p,hessIncrement);
h1 = -A0 \ Bq0cq0;                                             % -A\B(q0,conj(q0))
h2 = (2*1i*omega*speye(size(A0))-A0) \ Bq0q0;                  % (2iw-A)\B(q0,q0)

Cq0q0cq0 = multilinear3(q0,q0,conj(q0),x,p,ten3Increment);

Bq0h1    = multilinear2(q0, h1, x, p, hessIncrement);	       % +2B(h1,q0)
Bcq0h2   = multilinear2(conj(q0), h2, x, p, hessIncrement);	   % + B(conj(q0), h2)

coef = real(p0'*(Cq0q0cq0 + 2*Bq0h1 + Bcq0h2))/2.0;