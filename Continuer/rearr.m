function [x,p,k] = rearr(x0)
%
% [x,p] = rearr(x0)
%
% Rearranges x0 into coordinates (x) and parameters (p)
global cds
p = cds.P0;
p(cds.ActiveParams) = x0(cds.ncoo+(1:cds.nap));
x = x0(1:cds.ncoo);
k = x0(end);