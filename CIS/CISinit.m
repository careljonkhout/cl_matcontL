% [Q0,T0,evl,NSub,NUnstable] =  CISinit(A0, resizeable_flag, sparse_flag, ...
%                             NSub, NExtra, ...
%                             NUnstable, NStableRef, ...
%                             solver_opt, TwoParCont)
%
% Compute an orthonormal basis Q0 for the invariant subspace
% corresponding to the NSub least stable eigenvalues of A0.
%
% If resizeable_flag is true, then the size of the space NSub may
% be adjusted.  Then on output, NSub is the smallest number such that
%
%   NSub >= NSubRef = NUnstable + NStableRef
%
% and eigenvalues NSub and NSub + 1 (in descending order by real part)
% have real parts separated by a certain margin.  The code must compute
% at least (input NSub) + NExtra eigenvalues.
%
% If the input value for NUnstable is nonnegative and is different
% from the computed dimension of the unstable space, CISinit will exit
% with an error.  Also, if NSub is undefined or is greater than
% NSubRef + 2, then CISinit will exit with an error.

% JH: Modified second argument.  Possibly not needed. 10/12/06
% function [Q0,T0,evl,NSub,NUnstable] = CISinit(A0, resizeable_flag, sparse_flag, ...
%     NSub, NExtra, NUnstable, NStableRef, ...
%     solver_opt)
%

% JH: 06/07/07 There were problems with eigschurs_transform returning
% empty so several checks were added  to ensure that a subspace
% calculation would be possible

%function [Q0,T0,evl,NSub,NUnstable] = CISinit(A0, SparseSolvers, ...
 %   NSub, NExtra, NUnstable, NStableRef, solver_opt)
function [Q0, T0, evl, NSub, NUnstable] = CISinit(A0, NSub, NUnstable)
print_diag(5,'In CISinit\n');
%global cds cis contopts % MP
global cds contopts % MP

SparseSolvers  = contopts.CIS_SparseSolvers;
resizeable_flag= contopts.CIS_resizeable_flag; 
%NUnstable     = contopts.CIS_NUnstable;
MaxUnstable    = contopts.CIS_MaxUnstable;
NStableRef     = contopts.CIS_NStableRef;
NExtra         = contopts.CIS_NExtra;
%NUnstableGuess = contopts.CIS_NUnstableGuess; % DV MP

%                                                 % MP 2018 test 
if (NSub ~= NUnstable + NStableRef) & (NUnstable ~= -1)
   NSub 
   NUnstable
   NStableRef
   error('must have: NSub = NUnstable + NStableRef');
end
%                                                 % MP 2018 test 

%%JH: 5/14/2007 - Fixing potential error with unresolved overlap error in
%%first step.

NUnstableGuess = MaxUnstable; % MP
if NUnstable == -1 % calculate NUnstable for first step MP
    done = 0;
    while ~done

        %if NUnstableGuess > MaxUnstable || NUnstableGuess < 1
         %   %error('Cannot find sufficiently Small Unstable Space\n')
         %   error('must have: 1 < NUnstableGuess < MaxUnstable\n')   % MP
        %end

        if ~SparseSolvers
            [Q0,T0] = schur(full(A0));
            [~,Tt] = rsf2csf(Q0,T0);
            evlt    = diag(Tt);
        else
            if MaxUnstable >= length(A0)  %MP
                %error('Cannot find sufficiently Small Unstable Space\n')
                error('must have: MaxUnstable < length(A0)\n')   % MP
            end
             [Q,T]   = eigschurs_transform(A0, NUnstableGuess);
            [~,Tt] = rsf2csf(Q,T);
            evlt    = diag(Tt);
            if isempty(Q)
                NUnstableGuess = NUnstableGuess - 1;
                continue
            end
        end
        % Sort the eigenvalues in descending order by real part
        real_evl          = real(evlt);
        [real_evl, ~] = sort(-real_evl);
        real_evl          = -real_evl;
        %evl               = evlt(sorti);

        if isequal(cds.curve,@limitpointL) || isequal(cds.curve,@hopfL)
            [~,k] = min(abs(real_evl));
            nonzero_evl = real_evl;
            %NUnstable_t = sum(nonzero_evl > -1.0e5*eps);  % NASA
            NUnstable = k-1; %JH more accurate, ignores zero evals
        else
            NUnstable = sum(real_evl > -1.0e5*eps);  % NASA
        end

        %if NUnstable > NUnstableGuess - 4
         %   NUnstableGuess = NUnstableGuess + 5;
        %else
            done = done+1;
        %end
    end
    if resizeable_flag                       % MP 2018  
         NSub = NUnstable + NStableRef;% No initial NSub.  Smallest First guess
    end                                      % MP 2018       
    %cds.NUnstableGuess = NUnstable + 2; %% DV
    %cis.NUnstableGuess = NUnstable + 2; %% MP ??????????
end
%JH - End 5/14/2007 - NUnstable calculated for first step

% Do the initial (partial) eigendecomposition

done = 0;
Nevl = Inf;
Nevl = NSub + NExtra + 2;
while ~done
    if Nevl < NUnstable
        error('Cannot find sufficiently Small Unstable Space\n');
    end
    if ~SparseSolvers
        [Q0,T0] = schur(full(A0));
        [~,Tt] = rsf2csf(Q0,T0);
        evlt    = diag(Tt);
        if isempty(Q0) % DV: Can never happen? MP 7-2017 remove commrnts
            Nevl = Nevl-1; % MP 7-2017 remove commrnts
        else               % MP 7-2017 remove commrnts
             done = 1;     % MP 7-2017 remove commrnts
        end                % MP 7-2017 remove commrnts
    else
        [Q,T]   = eigschurs_transform(A0, Nevl);
        [~,Tt] = rsf2csf(Q,T);
        evlt    = diag(Tt);
        if isempty(Q)
            Nevl = Nevl-1;
        else
            done = 1;
        end
    end

end
% Sort the eigenvalues in descending order by real part
real_evl          = real(evlt);                       %% MF
[real_evl, sorti] = sort(-real_evl);                  %% MF
real_evl          = -real_evl;

if isequal(cds.curve,@limitpointL) || isequal(cds.curve,@hopfL)
    [~,k] = min(abs(real_evl));
    %NUnstable_t = k-1; % MP
    NUnstable = k-1; % MP
else
    %NUnstable_t = sum(real_evl > -1.0e5*eps);  % MP
    NUnstable = sum(real_evl > -1.0e5*eps);  % MP
end

%if NUnstable >= 0 && NUnstable ~= NUnstable_t
%    error('NUnstable = %d; NUnstable_t = %d', NUnstable, NUnstable_t);
%end

% Update NSub if needed
if resizeable_flag                          % MP 2018  
    NSubRef  = NUnstable + NStableRef;  %Test

    gaps     = real_evl(NSubRef:end-1) - real_evl(NSubRef+1:end);
    gaps_idx = find(gaps > 1.0e2*eps) + NSubRef - 1;

    if isempty(gaps_idx) || gaps_idx(1) > NSubRef + 2
        error('Unresolved overlap: three real parts of eigenvalues too close');
    else
        NSub = min(gaps_idx(1),length(sorti));
    end
    print_diag(2,'(Updated) NSub = %d NUnstable = %d\n',NSub, NUnstable);
end

% Reorder (partial) Schur form
if ~SparseSolvers
    % old [Q0, T0, WR, WI, s, sep] = mexdtrsen(Q0, T0, sorti(1:NSub)); %MF 2/21/2013
    E = ordeig(T0);  %MF 2/21/2013
    EE = E(sorti);  %MF 2/21/2013
    E_NSub = EE(NSub);  %MF 2/21/2013
    [Q0, T0] = ordschur(Q0, T0, E >= E_NSub);   %MF 2/21/2013    
    T = T0; % DV :quick fix
else
     % old [U,  T,  WR, WI, s, sep] = mexdtrsen(eye(length(T)), T, sorti(1:NSub)); %MF 2/21/2013  
    E = ordeig(T);  %MF 2/21/2013
    EE = E(sorti);  %MF 2/21/2013
    E_NSub = EE(NSub);  %MF 2/21/2013
    [U,T] = ordschur(eye(length(T)),T,E >= E_NSub);   %MF 2/21/2013 
   
    Qhat = Q*U;                   % Reorder Q corresponding to T
    Q0   = Qhat(:,1:NSub);
    T0   = T(1:NSub,1:NSub);
end
% old evl = complex(WR,WI);   %MF 2/21/2013
evl = ordeig(T);   %MF 2/21/2013