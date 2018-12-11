function [Q,T] = eigschurs_transform(A, N)

global contopts

transform    = contopts.CIS_Ric_Transform;
cayley_shift = contopts.CIS_Ric_Cayley_Shift;

if strcmp(transform, 'invert') % DV 2018
    sigma = 0;
elseif strcmp(transform, 'cayley')
    sigma = 1;
    if isnumeric(cayley_shift)   % DV MP
        sigma = cayley_shift;
    end
end
    
    % use specified options
    opts.blsz  = contopts.CIS_Ric_schur_blsz;
    opts.nbls  = contopts.CIS_Ric_schur_nbls;
    opts.maxit = contopts.CIS_Ric_schur_maxit;
    opts.tol   = contopts.CIS_Ric_schur_tol;
    
opts.dispr = 0;
opts.k     = N; 

if strcmpi(transform, 'invert')
    % Shift and invert transformation
    opts.sigma     = 'SR';          % DV: is this correct?
    
    Alen = length(A);
    Ashift  = A - sigma*speye(Alen);
    [L,U,P,QQ] = lu(Ashift);
    Afun =@(X, L, U, P, QQ) QQ*(U\(L\(P*X)));
    
    [Q, ~, FLAG] = ahbschur(Afun, Alen, opts, L, U, P, QQ);    
    
elseif strcmpi(transform, 'cayley')
    % Cayley transformation
    opts.sigma = 'LR';              % DV: is this correct?
    
    Alen      = length(A);
    A1        = sigma*speye(Alen) - A;
    A2        = sigma*speye(Alen) + A;
    [L1,U1,P1,QQ1] = lu(sparse(A1));
    Afun =@(X, L, U, P, QQ, A2) QQ*(U\(L\(P*A2*X)));
    
    [Q, ~, FLAG] = ahbschur(Afun, Alen, opts, L1, U1, P1, QQ1, A2);    
    
else
    % No transformation
    opts.sigma = 'LR';
    
    [Q,~,FLAG] = ahbschur(A, opts);    
    
end

% check Schur transform and remove possibly applied transformation
That      = Q'*A*Q;
failed = check_schur(That);

if FLAG(1) || failed
    Q = []; T = []; return
end

[Qhat, T] = schur(That);
Q         = Q*Qhat;

% -- Check to make sure that this really looks like a partial Schur form
function failed = check_schur(T)
failed = 1;
if norm(tril(T,-2)) > 1e-4*norm(T)
    print_diag(4, ['Bad clustering behavior in eigschurs; ', ...
        'consider changing shift parameters']);
else
    failed = 0;
end