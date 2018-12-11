function [a3, b3] = nf_BT3_L(CISdata,x,p, varargin)
%
% compute normal form coefficients a3 and b3 for Bogdanov-Takens.

global contopts

A0 = CISdata.A;
[q0, q1, p1, p0] = getBTvectors(CISdata, varargin{:});

hessIncrement = (contopts.Increment)^(3.0/4.0);
ten3Increment = (contopts.Increment)^(3.0/5.0);

Bq0q0 = multilinear2(q0,q0,x,p,hessIncrement);
Bq0q1 = multilinear2(q0,q1,x,p,hessIncrement);
Bq1q1 = multilinear2(q1,q1,x,p,hessIncrement);

a2  = p1'*Bq0q0/2.0 ;
b2  = p0'*Bq0q0 + p1'*Bq0q1;

Cq0q0q0 = multilinear3(q0, q0, q0, x, p, ten3Increment);
Cq0q0q1 = multilinear3(q0, q0, q1, x, p, ten3Increment);

Bord = [A0, p1; q0' 0];
h20 = Bord \ [(2*a2*q1 - Bq0q0); 0]; 
h20 = h20(1:end-1);
delta0 = 2*p0'*Bq0q1 + p1'*Bq1q1 - 2*p0'*h20;
h20 = h20 + delta0*q0;  % DV: The shift in delta0*q0 is strictly speaking not needed here, but is consistent with the computation of a4 and b4
h11 = Bord \ [(b2*q1 + h20 - Bq0q1); 0];
h11 = h11(1:end-1);

Bh20q0 = multilinear2(h20, q0, x, p, hessIncrement);
Bh11q0 = multilinear2(h11, q0, x, p, hessIncrement);
Bh20q1 = multilinear2(h20, q1, x, p, hessIncrement);

a3 = 1/6*p1'*Cq0q0q0 + 1/2*p1'*Bh20q0 ...
     - a2/2*p1'*Bq1q1;
 
b3 = 1/2*p1'*Cq0q0q1 + p1'*Bh11q0 + 1/2*p0'*Bh20q1 ...
     + 1/2*p0'*Cq0q0q0 + 3/2*p0'*Bh20q0 ...
     - b2/2*p1'*Bq1q1 +a2*p0'*Bq1q1 ...
     -5*a2*p0'*h11;