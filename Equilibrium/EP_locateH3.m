%% ------------------Hopf Locator----------------------------------------
function pout = locateH3(p1, p2)
pout = [];
x1 = p1.x; v1 = p1.v;
x2 = p2.x; v2 = p2.v;
NSub = p1.CISdata.NSub;

print_diag(5,'In equilibriumL/locateH3\n');

global cds contopts

%JH: Bisection 9/1/06----------
%v = v1 + v2;
%v = v/norm(v);
normdx = inf;
normf  = inf;
%NBisec = contopts.NBisectionIters;
X1 = x1; V1 = v1; X2 = x2; V2 = v2;
distX1X2 = norm(X1 - X2);
%------------------------------

ndim         = cds.ndim;
Incr         = contopts.Increment;
MaxCorrIters = contopts.contL_Loc_MaxCorrIters; 
VarTolerance = contopts.contL_Loc_VarTolerance; 
FunTolerance = contopts.contL_Loc_FunTolerance;  

JacMatrix = contjac(x1);
A2 = JacMatrix(1:ndim-1, 1:ndim-1);
eds2.A2 = A2; 
% DV: T1 is computed at currpoint and not at x1, so recompute
CISdata = contCIS_step(A2, p1.CISdata);
T1 = CISdata.T;
Q1 = CISdata.Q;
Q1_partial = Q1(:, 1:NSub);
A1 = A2;         

Asub   = T1(1:NSub,1:NSub);  
[V,D1] = eig(Asub);
d      = diag(D1);

% find a complex eigenvalue with smallest absolute real part and set omega
% = its imaginary part
I     = find(imag(d) ~= 0);
if isempty(I)
    print_diag(1, 'Locate Hopf: All eigenvalues in start point are real \n');
    x = []; v = [];
    return
end
Rmin  = min(abs(real(d(I))));
I1    = find(abs(real(d)) == Rmin);
lamda = real(d(I1(1)));
omega = imag(d(I1(1)));				% get omega of Asub

kapa  = lamda*lamda + omega*omega;

[W,D2]= eig(Asub');
d     = diag(D2);
I     = find(imag(d) ~= 0);
Rmin  = min(abs(real(d(I))));
I2    = find(abs(real(d)) == Rmin);

% Input
% make the bordered system nonsigular
[Q,R,E] = qr([real(V(:,I1(1))) imag(V(:,I1(1)))]);
borders.v =  [real(V(:,I1(1))) imag(V(:,I1(1)))];

[Q,R,E] = qr([real(W(:,I2(1))) imag(W(:,I2(1)))]);
borders.w =  [real(W(:,I2(1))) imag(W(:,I2(1)))];

x = [x1;kapa];   % the variables are (u,alpha,kapa)

%data structure for constructing minimally augmented system dynamically
two_index = zeros(6,2);

two_index(1,:) = [1 2];
two_index(2,:) = [1 3];
two_index(3,:) = [1 4];
two_index(4,:) = [2 3];
two_index(5,:) = [2 4];
two_index(6,:) = [3 4];

transferForm = [1,1;1,2;2,1;2,2];

for i = 1:MaxCorrIters
    
    RED = Asub*Asub+kapa*eye(NSub);
    %RED
    
    %Step 1
    Bord  = [RED borders.w; borders.v' zeros(2)];
    bunit = [zeros(NSub,2);eye(2)];
    
    %Step 2
    vext = Bord\bunit;
    v1   = vext(1:NSub,1);
    v2   = vext(1:NSub,2);
    
    %Step 3
    wext = Bord'\bunit;
    w1   = wext(1:NSub,1);
    w2   = wext(1:NSub,2);
    
    % DV: Compute true left eigenvector of A^2 + kapa*I
    RED2 = A1*A1+kapa*speye(length(A1));
    Bord2 = [RED2' Q1_partial*borders.v; borders.w'*Q1_partial' zeros(2)];
    bunit2 = [zeros(length(A1), 2); eye(2)];
    wext2 = bordCIS1(Bord2, bunit2, 2);
    Qw1   = wext2(1:end-2,1);
    Qw2   = wext2(1:end-2,2);
    
    %Step 4
    borders.v = [vext(1:NSub,1)/norm(vext(1:NSub,1)),vext(1:NSub,2)/norm(vext(1:NSub,2))];
    borders.w = [wext(1:NSub,1)/norm(wext(1:NSub,1)),wext(1:NSub,2)/norm(wext(1:NSub,2))];
    
    %Step 5    
    Qv1  = Q1_partial * v1;
    Qv2  = Q1_partial * v2;
    QCv1 = Q1_partial * Asub * v1;
    QCv2 = Q1_partial * Asub * v2;
    
    x1 = x(1:ndim);
    JacMatrix_x = JacMatrix(:,1:end-1);
    
    contCISjacDirDerSquare_x1_Q1_Asubv1 = DirDerJacobian(x1,Q1_partial*Asub*v1,JacMatrix_x);
    contCISjacDirDerSquare_x1_Q1_Asubv2 = DirDerJacobian(x1,Q1_partial*Asub*v2,JacMatrix_x);
    contCISjacDirDerSquare_x1_Q1_v1     = DirDerJacobian(x1,Q1_partial*v1     ,JacMatrix_x);
    contCISjacDirDerSquare_x1_Q1_v2     = DirDerJacobian(x1,Q1_partial*v2     ,JacMatrix_x);
    
    gx(1,:) = -Qw1'*contCISjacDirDerSquare_x1_Q1_Asubv1 - Qw1'*A1*contCISjacDirDerSquare_x1_Q1_v1;
    gx(2,:) = -Qw1'*contCISjacDirDerSquare_x1_Q1_Asubv2 - Qw1'*A1*contCISjacDirDerSquare_x1_Q1_v2;
    gx(3,:) = -Qw2'*contCISjacDirDerSquare_x1_Q1_Asubv1 - Qw2'*A1*contCISjacDirDerSquare_x1_Q1_v1;
    gx(4,:) = -Qw2'*contCISjacDirDerSquare_x1_Q1_Asubv2 - Qw2'*A1*contCISjacDirDerSquare_x1_Q1_v2;
    
    f_x1 = JacMatrix_x;
    x2   = x(1:ndim);
    x2(end) = x2(end) + Incr;
    f_x2 = contjac(x2);
    f_x2 = f_x2(1:end,1:ndim-1);
    
    gp(1,1) = -Qw1'*(f_x2*QCv1 - f_x1*QCv1)/(Incr) - w1'*Asub*Q1_partial'*(f_x2*Qv1 - f_x1*Qv1)/(Incr);
    gp(2,1) = -Qw1'*(f_x2*QCv2 - f_x1*QCv2)/(Incr) - w1'*Asub*Q1_partial'*(f_x2*Qv2 - f_x1*Qv2)/(Incr);
    gp(3,1) = -Qw2'*(f_x2*QCv1 - f_x1*QCv1)/(Incr) - w2'*Asub*Q1_partial'*(f_x2*Qv1 - f_x1*Qv1)/(Incr);
    gp(4,1) = -Qw2'*(f_x2*QCv2 - f_x1*QCv2)/(Incr) - w2'*Asub*Q1_partial'*(f_x2*Qv2 - f_x1*Qv2)/(Incr);
    
    gk(1,1) = -w1'*v1;
    gk(2,1) = -w1'*v2;
    gk(3,1) = -w2'*v1;
    gk(4,1) = -w2'*v2;
    
    %Step 6
    x1  = x(1:ndim);
    f_x = JacMatrix;
    
    if i == 1
        A = [f_x,zeros(ndim-1,1)];
        
        testMatrix = [A;gx gp gk];
        [Index1, Index2] = hopf_locator_construction(testMatrix, two_index);
    end;
    
    B   = [f_x,zeros(ndim-1,1); gx(Index1,:),gp(Index1,:),gk(Index1,:); gx(Index2,:),gp(Index2,:),gk(Index2,:)];
    f1  = [feval(cds.curve_func, x1); vext(NSub+transferForm(Index1,1),transferForm(Index1,2)); vext(NSub+transferForm(Index2,1),transferForm(Index2,2))];
    
    %Step 7
    dx = bordCIS1(B,f1,2);  % dx = B\f1
    
    x = x - dx;
    
    % DV: Seems we do not need this either
    %      if (norm(x(1:ndim) - X1) > distX1X2) & (norm(x(1:ndim) - X2) > distX1X2)
    %          print_diag(3,'Locate Hopf, Newton Produced Divergent Solution\n');
    %          x=[];        v=[];
    %          return
    %      end
    
    kapa = x(end);	% kapa is considered as variable now
    
    % user defined norm
    normdx_old = normdx;
    normf_old  = normf;
    normdx = norm(dx);
    normf  = norm(f1);
    
    print_diag(2,'Locate Hopf: Newton Step: norm(dx)=%6.9f norm(f1)=%6.9f\n', normdx, normf);
    
    if normdx < VarTolerance && normf < FunTolerance
        % End Locator Function, accept point.
        break;
    end
    
    %JH: Newton converging slowly, exit and reduce stepsize
    %DV: Exit and reduce stepsize does note work like that anymore
    %DV: Try just once and do not worry about monotonous convergence
    %     if (normf > 0.9*normf_old) || isnan(normf)
    %         print_diag(3,'Locate Hopf: Residual Converging Too Slowly\n');
    %         x=[];        v=[];
    %         return
    %     end
    
    if i == MaxCorrIters
        %when fail to locate Hopf point, report message and return(FQ, 08/01/2006)
        x = [];     v = [];
        print_diag(3,'locateHopf: Failed to find x in %d Newton iterations\n',MaxCorrIters);
        return;
    end
    
    %Step 8 Update f_x, Asub
    x_temp  = x(1:ndim);
    f_x     = contjac(x_temp);
    A2      = f_x(:,1:ndim-1);
    eds2.A2 = A2;                              
    JacMatrix = f_x;
    [CISdata] = contCIS_step(A2, CISdata);
    if isempty(CISdata)
        print_diag(3,'Locate Hopf: contCIS_step failed to find Asub');
    else
        A1   = A2;
        Q1   = CISdata.Q;
        Q1_partial = Q1(:, 1:CISdata.NSub);
        T1   = CISdata.T;            
        Asub = T1(1:CISdata.NSub,1:CISdata.NSub);
    end
end

kapa = x(end);
x = x(1:ndim);
v = V1+V2;
v = v/norm(v);

pout.x = x;
pout.v = v;
pout.tvals = [];
pout.uvals = [];

pout.s.borders.v = borders.v;
pout.s.borders.w = borders.w;
pout.s.kapa = kapa;