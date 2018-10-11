function coef = nf_HH_L(CISdata,x,p)
%
% compute normal form coefficients for double-hopf.
%
NSub   = CISdata.NSub;
A0     = CISdata.A;
T0     = CISdata.T(1:NSub,1:NSub);
Q0     = pout.Q(:,1:NSub);

[V, D]  = eig(T0); 
[W, D2] = eig(T0'); 

d       = diag(D);
d2      = diag(D2);

ind = find(imag(d) > 0);
d = d(ind);
V = V(:, ind); 
[~, ind] = sort(real(d)); 
ev1 = d(ind(1));
ev2 = d(ind(2));
v1 = V(:, ind(1)); 
v2 = V(:, ind(2));
q0 = Q0 * v1;
q1 = Q0 * v2;
borders1.v =  [real(v1) imag(v1)];
borders2.v =  [real(v2) imag(v2)];

ind = find(imag(d2) > 0);
d2 = d2(ind);
W = W(:, ind); 
[~, ind] = sort(real(d2)); 
w1 = W(:, ind(1)); 
w2 = W(:, ind(2)); 
borders1.w =  [real(w1) imag(w1)];
borders2.w =  [real(w2) imag(w2)];

Bord1 = [A0' Q0*borders1.v; borders1.w'*Q0' 0];
bunit2 = [zeros(length(A0), 2); eye(2)];
wext = bordCIS1(Bord1, bunit2, 2);  
w1   = wext(1:cds.ncoo,1);
w2   = wext(1:cds.ncoo,2);
alpha=w1'*A0'*w2+1i*omega*(w1'*w2);
beta=-w1'*A0'*w1-1i*omega*(w1'*w1);
p0=alpha*w1+beta*w2;

Bord2 = [A0' Q0*borders2.v; borders2.w'*Q0' 0];
bunit2 = [zeros(length(A0), 2); eye(2)];
wext = bordCIS1(Bord2, bunit2, 2);  
w1   = wext(1:cds.ncoo,1);
w2   = wext(1:cds.ncoo,2);
alpha=w1'*A0'*w2+1i*omega*(w1'*w2);
beta=-w1'*A0'*w1-1i*omega*(w1'*w1);
p1=alpha*w1+beta*w2;

% find adjoint vectors
q0=q0/norm(q0);
p0=p0/(q0'*p0);
q1=q1/norm(q1);
p1=p1/(q1'*p1);

%% NF computation

hessIncrement = (contopts.Increment)^(3.0/4.0);
ten3Increment = (contopts.Increment)^(3.0/5.0);
%2nd order vectors
h1100 = -A0\multilinear2(q0,conj(q0),x,p,hessIncrement);                           % -A\B(q0,conj(q0))
h2000 = (2*ev0*speye(size(A0))-A0)\multilinear2(q0,q0,x,p,hessIncrement);		   % (2iw_0-A)\B(q0,q0)			        
h1010 = ((ev0+ev1)*speye(size(A0))-A0)\multilinear2(q0,q1,x,p,hessIncrement);	   % (i(w_0+w_1)-A)\B(q0,q1)
h1001 = ((ev0-ev1)*speye(size(A0))-A0)\multilinear2(q0,conj(q1),x,p,hessIncrement);% (i(w_0-w_1)-A)\B(q0,conj(q1))
h0020 = (2*ev1*speye(size(A0))-A0)\multilinear2(q1,q1,x,p,hessIncrement);		   % (2iw_1-A)\B(q1,q1)
h0011 = -A0\multilinear2(q1,conj(q1),x,p,hessIncrement);                           % -A\B(q1,conj(q1))

%3rd order vectors
G2100 = p0' * multilinear3(q0,q0,conj(q0),x,p,ten3Increment) / 2.0 ...      %  C(q0,q0,conj(q0))
      + p0' * multilinear2(q0,h1100,x,p,hessIncrement) ...                  %+2B(h1100,q0)
      + p0' * multilinear2(h2000,conj(q0),x,p,hessIncrement) / 2.0;         %+ B(h2000,conj(q0))
  
G1011 = p0' * multilinear3(q0,q1,conj(q1),x,p,ten3Increment) ...            %  C(q0,q1,conj(q1))
      + p0' * multilinear2(h0011,q0,x,p,hessIncrement) ...                  %+ B(q0,h0011)
      + p0' * multilinear2(h1001,q1,x,p,hessIncrement) ...                  %+ B(q1,h1001)
      + p0' * multilinear2(h1010,conj(q1),x,p,hessIncrement);               %+ B(conj(q1),h1010)
  
H1110 = p1' * multilinear3(q1,q0,conj(q0),x,p,ten3Increment) ...            %  C(q1,q0,conj(q0))
      + p1' * multilinear2(h1100,q1,x,p,hessIncrement) ...                  %+ B(q1,h1100)
      + p1' * multilinear2(conj(h1001),q0,x,p,hessIncrement) ...            %+ B(q0,conj(h1001))
      + p1' * multilinear2(h1010,conj(q0),x,p,hessIncrement);               %+ B(conj(q0),h1010)
  
H0021 = p1' * multilinear3(q1,q1,conj(q1),x,p,ten3Increment) / 2. ...       %  C(q1,q1,conj(q1))
      + p1' * multilinear2(q1,h0011,x,p,hessIncrement) ...                  %+2B(h0011,q1)
      + p1' * multilinear2(h0020,conj(q1),x,p,hessIncrement) / 2.0;         %+ B(h0020,conj(q1))
  
%coefficients
coef = [ G2100 G1011; H1110 H0021 ];