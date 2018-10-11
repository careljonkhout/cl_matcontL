function CISdata = contCIS_step(A2, CISdata1)  %DV 2018

global cds contopts  %MP

% --- Read out all cds options
h_min = cds.h_min;
h     = cds.h;

% Read out CISdata of 
A1  = CISdata1.A;
Q1  = CISdata1.Q;
T1  = CISdata1.T;
NSub = CISdata1.NSub;

% -- Compute candidate step
if ~contopts.CIS_SparseSolvers
    [Q2, T2, evl2, newton_steps] = CISstep(Q1,T1,A2,NSub);
else
    % echeck                             MP
    %if size(Q1,2) ~= NSub
    %  print_diag(3,'contCIS_step: size(Q1,2) ~= NSub\n');
    %  size(Q1,2)
    % NSub
    %pause
    % error('size(Q1,2) ~= NSub\n');
    %end                                   MP
    [Q2, T2, evl2, newton_steps] = CISstepSp(Q1,A2,A1);  % MP
    
end

if isempty(Q2)
    if h <= h_min
        cds.tfUpdate = 1;
    end
    CISdata = [];
    return
end

% Sort eigenvalues in descending order; get "block" eigenvalues
[evl2_r, evl2_l] = SortEvl(evl2, NSub);

% Detect overlaps (if NSub is less than the size of Jaacobian matrix) MP
unresolved  = 0;                                         % MP 2018
overlap     = 0;                                         % MP 2018
needs_adapt = 0;                                         % MP 2018
if contopts.CIS_DetectOverlap                % MP
    is_conj_r    = (imag(evl2_r(end)) ~= 0);
    is_conj_l    = (imag(evl2_l(1))   ~= 0);
    
    rleft1       =  real(evl2_l(1));
    rleft2       =  real(evl2_l(2 + is_conj_l));
    overlap      = (real(evl2_r(end)) < rleft1);
    overlap2     = (real(evl2_r(end)) < rleft2);
    
    % Tell user about overlaps
    if (overlap2),   print_diag(3,'nongeneric ');    end
    if (overlap ),   print_diag(3,'overlap: \n');    end
    
    unresolved  = overlap2;
    needs_adapt = (contopts.CIS_AdaptOnOverlap & overlap);
end                                              % MP

% -- Step acceptance and adaptation

if unresolved
    print_diag(3,'contCIS_step: Decreasing step size\n');
    CISdata = [];
    return
else            % accept candidate step, adapting subspace if needed
    % DV: create CIS data struct
    CISdata.A = A2;
    CISdata.Q = Q2;
    CISdata.T = T2;
    CISdata.evl_r = evl2_r;
    CISdata.evl_l = evl2_l;
    CISdata.NSub = NSub;
    CISdata.NUnstable = CISdata1.NUnstable;
    CISdata.overlap = overlap;
end