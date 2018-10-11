% 
% [Y,evl,newton_steps] = RicSolve(T0, Th, Nsub, solver_opt)
%
% Given T0 and Th, set up and then solve for Y
% Riccati equation: T22h*Y - Y*T11h = -E21 + Y*T12h*Y
%
% Inputs:
%  T0:        an n-by-n (in dense case) old block Schur
%       (with nsub-dimensional invariant subspace near the imaginary axis)
%  Th:  1) dense case: Th = Q(s)'*A(s+h)*Q(s), an n-by-n
%  Th:  2) sparse case: a block upper Hessenberg matrix
% Outputs:
%  Y:        an Nsub-by-(n-Nsub) matrix

function [Y,evl,i] = RicSolve(T0,Th,NSub)

global contopts
Euler          = contopts.CIS_Ric_Euler;
FunTol         = contopts.CIS_Ric_FunTolerance;
VarTol         = contopts.CIS_Ric_VarTolerance;
MaxNewtonIters = contopts.CIS_Ric_MaxNewtonIters;

%  Set up the coefficient matrices T22h, T11h, E21, T12h
%
T11h  = Th(1:NSub,     1:NSub);
T12h  = Th(1:NSub,     NSub+1:end);
E21   = Th(NSub+1:end, 1:NSub);
T22h  = Th(NSub+1:end, NSub+1:end);

T0_11 = T0(1:NSub,     1:NSub);
T0_22 = T0(NSub+1:end, NSub+1:end);

% Euler predict Y (solve Sylvester equation T0_22*Y - Y*T0_11 = -E21)
%
n = length(Th);
Y = zeros(n-NSub,NSub);
if Euler == 1
  [Y,evl] = Sylvsol(T0_22, T0_11, -E21);
  if Y == []
    print_diag(3,'RicSolve: Could not solve Euler predictor equation\n')
    return %JH 09/08/2006
  end
end

F = T22h*Y - Y*T11h + E21 - Y*T12h*Y;

%  Newton corrections
%
%for i = 1:solver_opt.MaxNewtonIters
for i = 1:MaxNewtonIters    

  % Do one Newton iteration
  [dY,evl] = Sylvsol(T22h - Y*T12h, T11h + T12h*Y, -F);

  if isempty(dY)
    print_diag(3,'RicSolve: Could not solve Newton corrector equation\n')
    y = []; %JH 09/08/2006
    return  %JH 09/08/2006
  end

  % Add the Newton increment to update matrix Y :
  Y = Y + dY;

  % Check whether errors reached user-supplied tolerances :
  F     = T22h*Y - Y*T11h + E21 - Y*T12h*Y;
  normY = norm(Y, 'fro');
  erdY  = norm(dY,'fro')/(1 + normY);    %??? check!!
  normF = norm(F, 'fro');

  %if normF < solver_opt.FunTol & erdY < solver_opt.VarTol
  if normF < FunTol & erdY < VarTol
     return;
  end

end
print_diag(3,'RicSolve: Converging too slowly\n')
Y = [];
