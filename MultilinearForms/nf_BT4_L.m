function [a4, b4] = nf_BT4_L(CISdata,x,p, varargin)
%
% compute normal form coefficients a3 and b3 for Bogdanov-Takens.

global contopts

A0 = CISdata.A;
[q0, q1, p1, p0] = getBTvectors(CISdata, varargin{:});

hessIncrement = (contopts.Increment)^(3.0/4.0);
ten3Increment = (contopts.Increment)^(3.0/5.0);
ten4Increment = (contopts.Increment)^(3.0/6.0);

Bq0q0 = multilinear2(q0,q0,x,p,hessIncrement);
Bq0q1 = multilinear2(q0,q1,x,p,hessIncrement);
Bq1q1 = multilinear2(q1,q1,x,p,hessIncrement);

a2  = p1'*Bq0q0/2.0 ;
b2  = p0'*Bq0q0 + p1'*Bq0q1;

Cq0q0q0 = multilinear3(q0, q0, q0, x, p, ten3Increment);
Cq0q0q1 = multilinear3(q0, q0, q1, x, p, ten3Increment);
Cq0q1q1 = multilinear3(q0, q1, q1, x, p, ten3Increment);

Bord = [A0, p1; q0' 0];

h20 = Bord \ [(2*a2*q1 - Bq0q0); 0];
h20 = h20(1:end-1);
delta0 = p0'*Bq0q1 + 1/2*p1'*Bq1q1 - p0'*h20;
h20 = h20 + delta0*q0; 
h11 = Bord \ [(b2*q1 + h20 - Bq0q1); 0];
h11 = h11(1:end-1);
h02 = Bord \ [(2*h11 - Bq1q1); 0];
h02 = h02(1:end-1);

Bh20q0 = multilinear2(h20, q0, x, p, hessIncrement);
Bh11q0 = multilinear2(h11, q0, x, p, hessIncrement);
Bh20q1 = multilinear2(h20, q1, x, p, hessIncrement);
Bh02q0 = multilinear2(h02, q0, x, p, hessIncrement);
Bh11q1 = multilinear2(h11, q1, x, p, hessIncrement);

a3 = 1/6*p1'*Cq0q0q0 + 1/2*p1'*Bh20q0 ...
    - a2/2*p1'*Bq1q1;

b3 = 1/2*p1'*Cq0q0q1 + p1'*Bh11q0 + 1/2*p1'*Bh20q1 ...
    + 1/2*p0'*Cq0q0q0 + 3/2*p0'*Bh20q0 ...
    - b2/2*p1'*Bq1q1 +a2*p0'*Bq1q1 ...
    -5*a2*p0'*h11;

Dq0q0q0q0 = multilinear4(q0, q0, q0, q0, x, p, ten4Increment);
Dq0q0q0q1 = multilinear4(q0, q0, q0, q1, x, p, ten4Increment);

h30 = Bord \ [(6*a3*q1 + 6*a2*h11 - 3*Bh20q0 - Cq0q0q0); 0];
h30 = h30(1:end-1);
gamma0 = p0'*(h30 + 2*a2*h02 + 2*b2*h11 - 2*Bh11q0 - Bh20q1 - Cq0q0q1);
gamma1 = p1'*(2*b2*h02 - Bh02q0 - 2*Bh11q1 - Cq0q1q1);  
delta1 = -gamma0 - gamma1;
h30 = h30 + delta1*q0;
h21 = Bord \ [(h30 + 2*b3*q1 + 2*a2*h02 + 2*b2*h11 - 2*Bh11q0 - Bh20q1 - Cq0q0q1); 0];
h21 = h21(1:end-1);
% delta2 =  % DV: we do not need h03 so we can omit this transformation
% h21 = h21 + delta2*q0;


h12 = Bord \ [(2*h21 + 2*b2*h02 - Bh02q0 - 2*Bh11q1 - Cq0q1q1); 0];
h12 = h12(1:end-1);
% h03 = Bord \ [(3*h12 - 3*Bh02q1 - Cq1q1q1); 0];
% h03 = h03(1:end-1);

Bh30q0  = multilinear2(h30,q0 ,x,p,hessIncrement);
Bh20h20 = multilinear2(h20,h20,x,p,hessIncrement);
Bh21q0  = multilinear2(h21,q0 ,x,p,hessIncrement);
Bh11h20 = multilinear2(h11,h20,x,p,hessIncrement);
Bh30q1  = multilinear2(h30,q1 ,x,p,hessIncrement);

Ch20q0q0 = multilinear3(h20, q0, q0, x, p, ten3Increment);
Ch20q0q1 = multilinear3(h20, q0, q1, x, p, ten3Increment);
Ch11q0q0 = multilinear3(h11, q0, q0, x, p, ten3Increment);

a4 = 1/24*(p1'*Dq0q0q0q0 + 6*p1'*Ch20q0q0 ...
     + 4*p1'*Bh30q0 + 3*p1'*Bh20h20) ...
     -1/2*a2*p1'*h21 - a3*p1'*h11;
 
h40 = Bord \ [(24*a4*q1 + 12*a2*h21 + 24*a3*h11 - 4*Bh30q0 - 3*Bh20h20 - 6*Ch20q0q0 - Dq0q0q0q0); 0];
h40 = h40(1:end-1);

b4 = 1/6*p1'*Dq0q0q0q1 ...
    + 1/2*p1'*Ch20q0q1 + 1/2*p1'*Ch11q0q0 ...
    + 1/2*p1'*Bh21q0 + 1/2*p1'*Bh11h20 ...
    + 1/6*p1'*Bh30q1 - 1/6*p1'*h40 - 1/2*b2*p1'*h21 ...
    - a2*p1'*h12 - a3*p1'*h02 - b3*p1'*h11;