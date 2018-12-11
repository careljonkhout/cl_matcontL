function coef = nf_BT_L(CISdata,x,p, varargin)
%
% compute normal form coefficients for Bogdanov-Takens.

% it is optional to provide vectors q0, q1, p1, p0 in varargin{}
% when these are not provided these are computed using CIS
global contopts

[q0, q1, p1, p0] = getBTvectors(CISdata, varargin{:});

hessIncrement = (contopts.Increment)^(3.0/4.0);
Bq0q0 = multilinear2(q0,q0,x,p,hessIncrement);
Bq0q1 = multilinear2(q0,q1,x,p,hessIncrement);

a2 = p1'*Bq0q0/2.0;
b2 = p0'*Bq0q0 + p1'*Bq0q1;
coef = real([a2 b2]);