% Example 4
%   x1'=x1 -x2
%   x2'=    x2 -x3
%   x3'=      t*x3
%   x1(1)=x2(1)=x3(1)=1
close all;clear classes;
opt.RelTol=1e-6;opt.AbsTol=1e-6;
% Jacobian is banded (Jacobimatrix Bandstruktur)
% Jacobian is given (Jacobimatrix wird explizit angegeben)
opt.JacobianLowerBandwidth=0;
opt.JacobianUpperBandwidth=1;
opt.Jacobian=@jac;
% Jacobiauswertung ist so billig, deshalb
opt.RecomputeJACFactor=-1;

[tG,xG,stats]=seulexMex(@f,[1,2],[1;1;1],opt);

t=linspace(1,2,100);
h1=exp((t-1)/2.*(1+t));
plot(t,h1,'g',tG,xG(:,1),'rx',tG,xG(:,2),'bx',tG,xG(:,3),'mx');

