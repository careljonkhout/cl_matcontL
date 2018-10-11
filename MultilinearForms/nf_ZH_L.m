function coef = nf_ZH_L(CISdata, x, p) %MP
%
% coef= nf_cp(odefile,jacobian,hessians,x,p,vext,wext,nphase)
% compute normal form coefficients for zero-hopf.
% Have not checked the symbolic part at all.

global contopts cds

% find null vectors
NSub   = CISdata.NSub;
A0     = CISdata.A;
T0     = CISdata.T(1:NSub,1:NSub);
Q0     = CISdata.Q(:,1:NSub);

[V, D]  = eig(T0); 
[W, D2] = eig(T0'); 

d       = diag(D);
d2      = diag(D2);

%% Find eigenvectors
% fold part
[~, ind] = min(abs(d));
v0 = V(:, ind);
q0 = Q0 * v0;

[~, ind] = min(abs(d2));
w0 = W(:, ind);

Bord0 = [A0' Q0*v0; w0'*Q0' 0];
bunit = [zeros(length(A0),1); 1];
p0 = Bord0 \ bunit;
p0 = p0(1:cds.ncoo);

% hopf part
ind = find(sign(imag(d))==1);
d = d(ind);
V = V(:, ind);
[~, ind] = min(real(d));
v1 = V(:, ind);
q1 = Q0 * v1;
borders.v =  [real(v1) imag(v1)];
ev1 = d(ind);

ind = find(sign(imag(d2))==1);
d2 = d2(ind);
W = W(:, ind);
[~, ind] = min(real(d2));
w1 = W(:, ind);
borders.w =  [real(w1) imag(w1)];

Bord2 = [A0' Q0*borders.v; borders.w'*Q0' zeros(2)];
bunit2 = [zeros(length(A0), 2); eye(2)];

wext = bordCIS1(Bord2, bunit2, 2);  
w1   = wext(1:cds.ncoo,1);
w2   = wext(1:cds.ncoo,2);

omega = imag(d2(ind));
alpha=w1'*A0'*w2+1i*omega*(w1'*w2);
beta=-w1'*A0'*w1-1i*omega*(w1'*w1);
p1=alpha*w1+beta*w2;

% normalize
q0 = q0 / norm(q0);
p0 = p0/(q0'*p0);
q1 = q1 / norm(q1);
p1 = p1/(q1'*p1);


%% NF computations
hessIncrement = (contopts.Increment)^(3.0/4.0);
ten3Increment = (contopts.Increment)^(3.0/6.0);

h1 = multilinear2(q0, q0, x, p, hessIncrement);
h2 = multilinear2(q0, q1, x, p, hessIncrement);
h3 = multilinear2(q1, conj(q1), x, p, hessIncrement);

G200 = p0' * h1 / 2.0;   % 1/2 <p0, B(q0,q0)>
H110 = p1' * h2;	     %     <p1, B(q0,q1)>
G011 = p0' * h3;	     %     <p0, B(q1,bar(q1))>

h020 = (2*ev1*speye(size(A0)) - A0) \ multilinear2(q1, q1, x, p, hessIncrement);
h200 = [A0 q0 ; p0' 0] \ [2*G200*q0 - h1 ; 0];	        %-A^{INV}( B(q0,q0) - <p0,B(q0,q0)>q0 )
h200 = h200(1:end-1);
h011 = [A0 q0 ; p0' 0] \ [  G011*q0 - h3 ; 0];
h011 = h011(1:end-1);
h110 = [(ev1*speye(size(A0)) - A0) q1 ; p1' 0] \ [h2 - H110*q1 ; 0];
h110 = h110(1:end-1);

G300 = p0' * multilinear2(q0,h200, x, p, hessIncrement) / 6.0 ...
     + p0' * multilinear3(q0,q0,q0,x,p,ten3Increment) / 2.0; 

G111 = p0' * multilinear3(q0,q1,conj(q1),x,p,ten3Increment) ...
     + p0' * multilinear2(q0, conj(h110), x, p, hessIncrement) ...
     + p0' * multilinear2(conj(q0), h110, x, p, hessIncrement) ...
     + p0' * multilinear2(q0, h011, x, p, hessIncrement);

H210 = p1' * multilinear3(q0, q0, q1, x, p, ten3Increment) / 2.0 ...
     + p1' * multilinear2(q0, h110, x, p, hessIncrement) ...
     + p1' * multilinear2(q1, h200, x, p, hessIncrement) / 2.0;

H021 = p1' * multilinear3(q0, q1, conj(q1), x, p, ten3Increment) / 2.0 ...
     + p1' * multilinear2(q1, h011, x, p, hessIncrement) ...
     + p1' * multilinear2(conj(q1), h020, x, p, hessIncrement) / 2.0;

B0 = G200;
C0 = G011;
E0 = real(H210 + H110 * (real(H021)/G011 - 3*G300/(2*G200) + G111/(2*G011)) - H021*G200/G011);
    
coef = [B0, C0, E0];