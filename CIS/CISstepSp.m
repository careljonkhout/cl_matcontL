% 
% [V, T, evl, newton_steps] = CISstepSp(V0, A, A0, NExtra, solver_opt)
%
% Build a refined approximate invariant subspace V of A and
% corresponding Rayleigh block T.  Use ordinary Arnoldi from
% ARPACK to get the chosen basis.  Also returns NExtra
% eigenvalues evl from the complementary subspace, and the number
% of Newton iterations newton_steps.

%function [V, T, evl, newton_steps] = ...
%         CISstepSp(V0, A, A0, NExtra, solver_opt)
function [V, T, evl, newton_steps] = CISstepSp(V0, A, A0)
        
global contopts
Select = contopts.CIS_Ric_SubspaceSelect;

if strcmpi(Select, 'eig')
    [V, T, evl, newton_steps] = CISstepSp_eig(V0, A, A0);
else
    [V, T, evl, newton_steps] = CISstepSp_ric(V0, A, A0);
end
  
% ----

function [V, T, evl, newton_steps] = CISstepSp_eig(V0, A, A0)

global contopts

Nsub      = size(V0,2); % ?
opts.disp = 0;
opts.v0   = V0 * ones(Nsub1);
Nevl      = Nsub + contopts.CIS_NExtra;

[V,Aproj] = eigschurs_transform(A, Nevl);
if isempty(V)                                               % MF
  warning('CISstepSp_eig: eigschurs_transform failed to find V')
  return                                                    % MF
end                                                         % MF
[Qt,Tt]   = rsf2csf(V,Aproj);
lambda    = diag(Tt);
target    = eig(V0'*A*V0);

lambda_m  = lambda * ones(1, length(target));
target_m  = target * ones(1, length(lambda));

[dist_m, id] = sort(min(abs(target_m - lambda_m'))); 
E = ordeig(Aproj);  %MF 2/21/2013
EE = E(id);  %MF 2/21/2013
E_Nsub = EE(Nsub);  %MF 2/21/2013  
[Qreord,T] = ordschur(eye(length(Aproj)),Aproj,E >= E_Nsub);   %MF 2/21/2013 
     
V            = V*Qreord(:,1:Nsub);
T            = T(1:Nsub, 1:Nsub);

evl          = [lambda(id(1:Nsub)); lambda(id(Nsub+1:end))];
newton_steps = 5;


% ----

function [V, T, evl, newton_steps] = CISstepSp_ric(V0, A, A0)  % MP

% -- Get a basis for the least stable eigenspace
%    using eigschurs.

global contopts

Nsub      = size(V0,2);   %? MP 

Nevl = Nsub + contopts.CIS_NExtra;


[V,Aproj] = eigschurs_transform(A, Nevl);
if isempty(V)                                               % MF
  print_diag(4, 'CISstepSp_ric: eigschurs_transform failed to find V');
  V= []; T=[]; evl=[]; newton_steps=[];
  return                                                    % MF
end                                                         % MF

% -- To pick a subspace, solve a dense Riccati equation using
%    a Newton iteration.

% Change basis again so the leading columns
% of V are the best approx to V0 in our projection space.

[CU,CS,CV] = svd(V'*V0);
V1h        = CU(:,1:Nsub)*CV'; 
V2h        = CU(:,Nsub+1:end);
Vh         = [V1h, V2h];
V          = V * Vh;

% DSB -- Should probably not hardwire this tolerance
if (1-CS(Nsub,Nsub) > contopts.CIS_MaxAngle)
  V = []; T = []; evl=[]; newton_steps=[];
  print_diag(4, 'CISstepSp_ric: Angle to projection subspace is too big');
  return;
end

% Solve for Y the Riccati equation T22h*Y - Y*T11h = -T21h + Y*T12h*Y
% and get eigenvalues from subspace (evl_r) plus a few of eigenvalues
% from the complementary subspace (evl_compl), evl = [evl_r; evl_compl]
Th = Vh' *  Aproj * Vh;
T0 = V'  *  A0    * V;
[Y,evl,newton_steps] = RicSolve(T0,Th,Nsub);

if isempty(Y)
  
    V= []; T=[]; evl=[]; newton_steps=[];
  print_diag(4, 'CISstepSp_ric: RicSolve failed to find Y');
  return;
end

% -- Compute an orthonormal basis V for span([I; Y])
%    such that \|V - V0\| is minimal.

Vh         = [eye(Nsub); Y];
[CU,~,CV]  = svd(Vh,0);
Vh         = CU*CV';
T          = Vh'*Th*Vh;               % Corresponding Schur block
V          = V * Vh;