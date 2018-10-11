% function [x,v] = locateLP(x1, x2, v1, v2, A1, Q1, T1)
function pout = EP_locateLP(p1, p2)
print_diag(5,'In equilibriumL/locateLP\n');
pout = [];
x1 = p1.x; v1 = p1.v;
x2 = p2.x; v2 = p2.v;

global cds contopts

cds.s = [];

ndim = cds.ndim;
Incr = contopts.Increment;
MaxCorrIters = contopts.contL_Loc_MaxCorrIters;
VarTolerance = contopts.contL_Loc_VarTolerance;
FunTolerance = contopts.contL_Loc_FunTolerance;
SparseSolvers = contopts.CIS_SparseSolvers;

f_x = contjac(x1);
A0  = f_x(1:ndim-1, 1:ndim-1);
T0 = p1.CISdata.T;
Q1 = p1.CISdata.Q;
NSub = p1.CISdata.NSub;
CISdata = p1.CISdata;

Q1_partial = Q1(:,1:NSub);
Asub = T0(1:NSub,1:NSub);
d    = eig(Asub);

% find a real eigenvalue with smallest absolute value and call it Rmin
I    = find(imag(d) == 0);
Rmin = min(abs((d(I))));
I1   = find(abs((d)) == Rmin);
RED  = Asub - d(I1(1))*eye(NSub);

% Input
borders.v = null(RED);
borders.w = null(RED');

x = x1;

normdx = inf;
normf  = inf;

X1 = x1; V1 = v1; X2 = x2; V2 = v2;
distX1X2 = norm(X1 - X2);

%------------------------------
for i = 1:MaxCorrIters
    %Step 1
    Bord  = [Asub borders.w; borders.v' 0];
    bunit = [zeros(NSub,1);1];
    
    %Step 2
    vext  = Bord\bunit;
    
    %Step 3
    wext  = Bord'\bunit;
    Q1_partial = Q1(:,1:NSub);
    
    % DV: To compute g_x we need the left eigenvector, Q*w is not good
    
    % Bord2 = [A0' Q1*borders.v; borders.w'*Q1' 0]; orig %MP 4-2016
    Bord2 = [A0' Q1_partial*borders.v; borders.w'*Q1_partial' 0];  %MP 4-2016
    bunit2 = [zeros(length(A0), 1); 1];
    wext2 = bordCIS1(Bord2, bunit2, 1);
    Qw = wext2(1:end-1);
    
    %Step 4
    borders.v = vext(1:NSub)/norm(vext(1:NSub));
    borders.w = wext(1:NSub)/norm(wext(1:NSub));
    
    %Step 5
    Qv = Q1_partial * vext(1:NSub);
    
    x1 = x;
    x1(1:ndim-1) = x1(1:ndim-1) - Incr/norm(Qv)*Qv;
    f_x1 = contjac(x1);
    x2 = x;
    x2(1:ndim-1) = x2(1:ndim-1) + Incr/norm(Qv)*Qv;
    f_x2 = contjac(x2);
    
    g_x = -(Qw'*f_x2 - Qw'*f_x1) * norm(Qv)/(2*Incr); % DV: minus sign was missing
    
    %Step 6
    B   = [f_x; g_x];
    
    f1  = [feval(cds.curve_func, x); vext(end)]; %orig %MP 4-2016
    %    f1_0 = feval(cds.curve_func, x);              %MP 4-2016
    %    if ~SparseSolvers                              %MP 4-2016
    %       f1 = [f1_0 vext(end)];                     %MP 4-2016
    %       f1=f1';                                    %MP 4-2016
    %    elseif SparseSolvers                            %MP 4-2016
    %       f1 = [f1_0; vext(end)];                    %MP 4-2016
    %    end                                           %MP 4-2016
    
    %Step 7
    dx = bordCIS2(B,f1,1);  % dx = B\f1
    
    x = x - dx;
    
    if (norm(x(1:ndim) - X1) > distX1X2) & (norm(x(1:ndim) - X2) > distX1X2)
        print_diag(3,'Locate Hopf, Newton Produced Divergent Solution\n');
        x=[];        v=[];
        return
    end
    
    % user defined norm
    normdx_old = normdx;
    normf_old  = normf;
    normdx = norm(dx);
    normf  = norm(f1);
    
    print_diag(2,'Locate LP, Newton Step: norm(dx)=%6.9f norm(f1)=%6.9f\n', normdx, normf);
    if normdx < VarTolerance && normf < FunTolerance
        break;
    end
    
    %JH: Newton converging slowly, exit and reduce stepsize
    if ((normf > normf_old) || isnan(normf))
        print_diag(3,'Locate LP: Residual Converging Too Slowly\n');
        x=[];       v = [];
        return
    end
    if i == MaxCorrIters
        print_diag(3,'Locate LP: Failed to find x with %d Newton iterations\n',MaxCorrIters);
        x = [];     v = [];
        return;
    end
    
    %Step 8                 Update f_x, Asub
    f_x = contjac(x);
    A2 = f_x(:,1:ndim-1);
    CISdata = contCIS_step(A2, CISdata);
    if isempty(CISdata)                                      % CIS failed
        print_diag(3,'Locate LP: contCIS_step failed to find Asub\n');
    else
        T1   = CISdata.T;
        Q1   = CISdata.Q;
        Q1_partial = Q1(:,1:CISdata.NSub);
        Asub = T1(1:CISdata.NSub,1:CISdata.NSub);
    end
    
end

v = v1 + v2;
v = v/norm(v);

pout.x = x;
pout.v = v;
pout.R = [];
pout.tvals = [];
pout.uvals = [];