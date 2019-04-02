function df=jac(t,x)
% Jacobian
df={[-1,-1],[1,1,t]};
% also possible:
% df=[ NaN,-1,-1; 1,1,t];
