% 

% [Q,T,evl,newton_steps] = CISstep(Q0, T0, A, NSub, solver_opt)
%
% Continue the block Schur form (Q0,T0) at A0 to a block Schur form
% (Q,T) at A.  NSub is the dimension of the leading block.
% evl(1:NSub) are eigenvalues for the leading block of T; evl(Nsub:end)
% are eigenvalues for the trailing block.  newton_steps is a count
% of the number of Newton iterations to converge.

function [Q,T,evl,newton_steps] = CISstep(Q0,T0,A,NSub)

global contopts
Select = contopts.CIS_Ric_SubspaceSelect;        

if strcmpi(Select, 'eig') 
    [Q, T, evl, newton_steps] = CISstep_eig(Q0,T0,A,NSub);
else
    [Q, T, evl, newton_steps] = CISstep_ric(Q0,T0,A,NSub);
end

% ----

function [Q,T,evl,newton_steps] = CISstep_eig(Q0,T0,A,NSub)

global contopts                              

PartialQ = contopts.CIS_Ric_PartialQ;        

[Q, T ] = schur(full(A));
[Qt,Tt] = rsf2csf(Q,T);

lambda  = diag(Tt);
target  = eig(Q0(:,1:NSub)'*A*Q0(:,1:NSub));

lambda_m  = lambda * ones(1, length(target));
target_m  = target * ones(1, length(lambda));

[dist_m, id] = sort(min(abs(target_m - lambda_m')));
%% [Q, T]       = mexdtrsen(Q, T, id(1:NSub));  %MF 2/21/2013
E = ordeig(T);  %MF 2/21/2013
EE = E(id);  %MF 2/21/2013
E_NSub = EE(NSub);  %MF 2/21/2013
[Q, T] = ordschur(Q, T, E >= E_NSub);   %MF 2/21/2013  

evl          = [lambda(id(1:NSub)); lambda(id(NSub+1:end))];
newton_steps = 5;

if PartialQ
  Q = Q(:,1:NSub+1);
  T = T(1:NSub+1,1:NSub+1);
end


% ----

function [Q,T,evl,newton_steps] = CISstep_ric(Q0,T0,A,NSub)
% Solve for Y the Riccati equation: T22h*Y - Y*T11h = -E21 + Y*T12h*Y

Th                   = Q0'*A*Q0;
[Y,evl,newton_steps] = RicSolve(T0,Th,NSub);

if isempty(Y)    
  print_diag(4, 'CISstep_ric: RicSolve failed to find Y');
  Q = []; T = []; evl = []; newton_steps = 0; 
  return;
end

% Normalize the solution to find next point (Q,T) on the branch

Q01        = Q0(:,1:NSub);
Q02        = Q0(:,NSub+1:end);
[U1,~,R1] = svd(Q01 + Q02*Y , 0);
[U2,~,R2] = svd(Q02 - Q01*Y', 0);
Q          = [U1*R1', U2*R2'];
T          = Q'*A*Q;