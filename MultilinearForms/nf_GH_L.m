function [coef] = nf_GH_L(CISdata,x,p,kapa,borders)
%
%calculates the second lyapunov coefficient
% 

global cds contopts

NSub   = CISdata.NSub;
A0     = CISdata.A;
T0     = CISdata.T(1:NSub,1:NSub);
Q0     = CISdata.Q(:,1:NSub);

if nargin <= 3 % compute kapa if not given
    [V, D] = eig(Asub);
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
    
    [W, D] = eig(Asub');
    d = diag(D);
    I     = find(imag(d) ~= 0);
    Rmin  = min(abs(real(d(I))));
    I2    = find(abs(real(d)) == Rmin);
    borders.w =  [real(W(:,I2(1))) imag(W(:,I2(1)))];
    
    kapa = abs(d(I2(1)))^2;
end

% right eigenvector
NSub = size(T0, 1);
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

hessIncrement = (contopts.Increment)^(3.0/4.0);
ten3Increment = (contopts.Increment)^(3.0/5.0);
ten4Increment = (contopts.Increment)^(3.0/6.0);
ten5Increment = (contopts.Increment)^(3.0/7.0);

%2nd order vectors
h20 = (2i*omega*speye(size(A0))-A0)\multilinear2(q0,q0,x,p,hessIncrement);	% (2iw-A)\B(q0,q0)
h11 = -A0\multilinear2(q0,conj(q0),x,p,hessIncrement);			            % -A\B(q0,conj(q0))
%3rd order vectors
h30 = multilinear3(q0,q0,q0,x,p,ten3Increment);				                %  C(q0,q0,conj(q0))
h30 = h30 + 3*multilinear2(q0,h20,x,p,hessIncrement);			            %+3B(h20,q0)  
h30 = (3i*omega*speye(size(A0))-A0)\h30;
h21 = multilinear3(q0,q0,conj(q0),x,p,ten3Increment);			            %  C(q0,q0,conj(q0))
h21 = h21 + 2*multilinear2(q0,h11,x,p,hessIncrement);			            %+2B(h11,q0)
h21 = h21 + multilinear2(h20,conj(q0),x,p,hessIncrement);		            %+ B(h20,conj(q0))
g21 = p0'*h21/2.0;
h21 = [A0-1i*omega*speye(size(A0)) q0; p0' 0 ]\[ 2*g21*q0-h21 ; 0];
h21 = h21(1:end-1);
%4th order vectors
h31 = multilinear4(q0,q0,q0,conj(q0),x,p,ten4Increment);			        %  D(q0,q0,q0,conj(q0))
h31 = h31 + 3*multilinear3(q0,q0,h11,x,p,ten3Increment);			        %+3C(q0,q0,h11)
h31 = h31 + 3*multilinear3(q0,conj(q0),h20,x,p,ten3Increment);		        %+3C(q0,conj(q0),h20)
h31 = h31 + 3*multilinear2(h20,h11,x,p,hessIncrement);			            %+3B(h20,h11)
h31 = h31 + 3*multilinear2(h21,q0,x,p,hessIncrement);			            %+3B(h21,q0)
h31 = h31 +   multilinear2(h30,conj(q0),x,p,hessIncrement);		            %+ B(h30,conj(q0))
h31 = (2i*omega*speye(size(A0))-A0)\(h31 - 3*g21*h20);  
h22 = multilinear4(q0,q0,conj(q0),conj(q0),x,p,ten4Increment);		        %  D(q0,q0,conj(q0),conj(q0))
h22 = h22 + 4*multilinear3(q0,conj(q0),h11,x,p,ten3Increment);		        %+4C(q0,conj(q0),h11)
h22 = h22 + 2*real(multilinear3(conj(q0),conj(q0),h20,x,p,ten3Increment));	%+2*Re(C(q0,q0,h02))
h22 = h22 + 4*real(multilinear2(h21,conj(q0),x,p,hessIncrement));		    %+2*Re(B(h21,conj(q0)))
h22 = h22 + 2*multilinear2(h11,h11,x,p,hessIncrement);			            %+2B(h11,h11)
h22 = -A0\(h22 + multilinear2(h20,conj(h20),x,p,hessIncrement));	        %+ B(h20,h02)
%5th order rhs
h32 = multilinear5(q0,q0,q0,conj(q0),conj(q0),x,p,ten5Increment);	        %  E(q0,q0,q0,conj(q0),conj(q0))
h32 = h32 + 6*multilinear4(q0,q0,conj(q0),h11,x,p,ten4Increment);	        %+6D(q0,q0,conj(q0),h11)
h32 = h32 + 3*multilinear4(conj(q0),conj(q0),q0,h20,x,p,ten4Increment);	    %+3D(conj(q0),conj(q0),q0,h20)
h32 = h32 +   multilinear4(q0,q0,q0,conj(h20),x,p,ten4Increment);	        %+ D(q0,q0,q0,h02)
h32 = h32 + 6*multilinear3(h11,h11,q0,x,p,ten3Increment);		            %+6C(h11,h11,q0)
h32 = h32 + 6*multilinear3(conj(q0),h20,h11,x,p,ten3Increment);		        %+6C(conj(q0),h20,h11)
h32 = h32 + 6*multilinear3(conj(q0),q0,h21,x,p,ten3Increment);		        %+6C(conj(q0),q0,h21)
h32 = h32 + 3*multilinear3(q0,h20,conj(h20),x,p,ten3Increment);		        %+3C(conj(q0),h20,h02)
h32 = h32 + 3*multilinear3(q0,q0,conj(h21),x,p,ten3Increment);		        %+3C(q0,q0,h12)
h32 = h32 +   multilinear3(conj(q0),conj(q0),h30,x,p,ten3Increment);	    %+ C(conj(q0),conj(q0),h30)
h32 = h32 + 6*multilinear2(h21,h11,x,p,hessIncrement);			            %+6B(h21,h11)
h32 = h32 + 3*multilinear2(h22,q0,x,p,hessIncrement);			            %+3B(h22,q0)
h32 = h32 + 3*multilinear2(h20,conj(h21),x,p,hessIncrement);		        %+3B(h12,h20)
h32 = h32 + 2*multilinear2(h31,conj(q0),x,p,hessIncrement);		            %+2B(h31,conj(q0))
h32 = h32 +   multilinear2(h30,conj(h20),x,p,hessIncrement);		        %+ B(h30,h02)
coef = real(p0'*h32/12);