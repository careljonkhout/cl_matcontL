% function [x,v] = locateBP(id, x1, x2, v1, v2)
function pout = EP_locateBP(p1, p2)
print_diag(5,'In equilibriumL/locateBP\n');
%(MP) reliability: combine Newton with bisection

global cds contopts  % MP

pout = [];
x1 = p1.x; v1 = p1.v;
x2 = p2.x; v2 = p2.v;

x1init = x1;
x2init = x2;

NumEvals = p1.CISdata.NSub + contopts.CIS_NExtra; % MP
ndim   = cds.ndim;
Incr   = contopts.Increment;
BordEr = contopts.contL_EQ_SingLinSystTol;
MaxCorrIters = contopts.contL_Loc_MaxCorrIters;
VarTolerance = contopts.contL_Loc_VarTolerance;
FunTolerance = contopts.contL_Loc_FunTolerance;

normdx = inf;
normf  = inf;

x  = 0.5*(x1+x2);
vx = v1 + v2;
vx = vx/norm(vx);

J = contjac(x);
if ~issparse(J)
    [V,evl] = eig([J;vx']);
    [d,evl1]= eig(J(:,1:ndim-1)');
    evl_vect = diag(evl);
    index = find(abs(evl_vect) == min(abs(evl_vect)));
    v = V(:,index(1));
else
    
    shift = 1e-2;
    B = [J;vx'];
    lB = size(B);
    R = rand(lB(1),1);
    R = R/norm(R);
    B = B - shift*speye(lB);
    v = bordCIS2(B,R,1);% v = B\R
    v = v/norm(v);
    
    if v ~= real(v)
        v = real(v);
    end
    opt.disp = 0;
    [d,evl1] = eigs(J(:,1:ndim-1), min(6,ndim-1), 'SM', opt);  % dsb %Question here about "6" eigenvalues, why this magic number?
    
end

bord.v = vx/norm(vx);
bord.w = v;

[y2,j] = min(abs(diag(evl1)));
bord.d = d(:,j);
i   = 0;

f1 = [zeros(ndim-1,2); eye(2)];
f3 = [zeros(ndim,1); 1];
vp  = cell(2);

bord.v0 = zeros(size(bord.v));
bord.w0 = zeros(size(bord.w));
i_improve = 0;
firstnewtonstep = 1;
while i <= MaxCorrIters
    
    %Step 1
    bord.w  = bord.w - (bord.w'*bord.v)*bord.v;
    bord.w  = bord.w/norm(bord.w);
    bord.W  = [bord.v bord.w];
    er_w    = norm(bord.w - bord.w0);
    er_v    = norm(bord.v - bord.v0);
    bord.v0 = bord.v;
    bord.w0 = bord.w;
    
    B = [J       bord.d
        bord.W' zeros(2,1)];
    
    if (i_improve < 1) && (er_w > BordEr || er_v > BordEr)
        % improve bord.d, bord.W
        i_improve = i_improve + 1;
        %Step 2
        [vg,fail] = bordCIS2(B,f1,2);  % vg = B\f1
        if fail
            % Newton seems to fail, use bisection
            print_diag(3,'Locate BP: Newton Fail: bordCIS2(B,f1,2) failed.\n')
            break
        end
        vp{1} = vg(1:ndim,1);
        vp{2} = vg(1:ndim,2);
        
        %Step 4
        [psig,fail] = bordCIS2(B',f3,2);  % psig = B'\f3
        if fail
            print_diag(3,'Locate BP: Newton Fail: bordCIS2(B(trans),f3,2) failed.\n')
            break
        end
        psi  = psig(1:ndim-1);
        
        %Step 8
        bord.d = psi/norm(psi);
        
        %Step 9
        bord.v = vp{1}/norm(vp{1});
        bord.w = vp{2};
        
    else
        %Step 2, Step 3
        f2  = -[feval(cds.curve_func, x); zeros(2,1)];
        f12 = [f1 f2];
        [vgh,fail] = bordCIS2(B,f12,2);  % vgh = B\f12
        if fail
            print_diag(3,'Locate BP: Newton Fail: bordCIS2(B,f12,2) failed.\n')
            break
        end
        
        %JH trying to avoid degeneracies
        nulldim = size(vgh);
        if nulldim(2) ~= 3
            print_diag(3,'Locate BP: Null Dimension ~= 3, Degenerate Case\n')
            return
        end
        
        vp{1} = vgh(1:ndim,1);
        vp{2} = vgh(1:ndim,2);
        xh   = vgh(1:ndim,3);
        muh  = vgh(ndim+1,3);
        
        %Step 4
        %tic%JH Testing NF Times
        [psig,fail] = bordCIS2(B',f3,2);  % psig = B'\f3
        if fail
            print_diag(3,'Locate BP: Newton Fail: bordCIS2(B(trans),f3,2) failed.\n')
            break;
        end
        psi  = psig(1:ndim-1);
        g    = psig(ndim:ndim+1);
        
        %Step 5
        dv{1} = Incr*vp{1}/norm(vp{1});
        dv{2} = Incr*vp{2}/norm(vp{2});
        dxh   = Incr*xh/norm(xh);
        
        for j = 1:2
            if norm(xh,inf) < eps
                eta(j,1) = 0;
            else
                eta(j,1) = norm(vp{j}) * norm(xh) / (Incr*Incr) * psi' * ...
                    (feval(cds.curve_func, x+dv{j}+dxh) - ...
                    feval(cds.curve_func, x+dv{j}) -     ...
                    feval(cds.curve_func, x+dxh) +       ...
                    feval(cds.curve_func, x));
            end
            
            M(j,j) = norm(vp{j})^2/(Incr*Incr)*psi'*...
                (  feval(cds.curve_func, x+dv{j}) - ...
                2*feval(cds.curve_func, x      ) + ...
                feval(cds.curve_func, x-dv{j}));
            
        end
        % new
        %print_diag(0,'NF time: ')%JH Testing NF Times
        
        M(1,2) = norm(vp{1})*norm(vp{2})/(Incr*Incr)*psi'*...
            (feval(cds.curve_func, x+dv{1}+dv{2}) - ...
            feval(cds.curve_func, x+dv{1}) - ...
            feval(cds.curve_func, x+dv{2}) + ...
            feval(cds.curve_func, x));
        M(2,1) = M(1,2);
        %toc%JH Testing NF Times
        %Step 6
        ksi = M\(g - eta);
        
        %Step 7
        dx = xh + ksi(1)*vp{1} + ksi(2)*vp{2};
        
        %Step 8 (modified for line search -- DSB)
        bord.d = psi/norm(psi);
        
        normdx_old = normdx;
        normf_old  = normf;
        
        x = x + dx;
        
        %         if norm(x-x1) > cds.h || norm(x-x2) > cds.h
        %             print_diag(3,'Locate BP: Newton Produced Divergent Solution\n');
        %             x=[];  v=[];
        %             return
        %         end
        
        f = feval(cds.curve_func, x) + muh*bord.d;
        
        normf  = normU(f);
        normdx = normU(dx);
        
        print_diag(2,'Locate BP: Newton Step: norm(dx)=%6.9f norm(f1)=%6.9f\n', normdx, normf);
        
        %Original Stopping Condition
        
        %% normdx                                                %%%
        %%VarTolerance
        %%normf
        %%FunTolerance
        %%%pause                                              %%%
        
        
        if normdx < VarTolerance && normf < FunTolerance
            v = v1 + v2;
            v = v/norm(v);
            
            pout.x = x;
            pout.v = v;
            pout.tvals = [];
            pout.uvals = [];
            
            v1 = v;
            v2 = vp{2} - v1*(vp{2}'*v1);
            v2 = v2/norm(v2);
            
            switch contopts.contL_EQ_BranchingMethod
                case 0  %Normal Form Coefficients
                    pout.s.v1 = v1;
                    pout.s.v2 = v2 - 0.5*(M(2,2)/M(1,2))*v1;
                    pout.s.v2 = pout.s.v2/norm(pout.s.v2);
                case 1  %Perpendicular Guess
                    pout.s.v1 = v1;
                    pout.s.v2 = v2;
                case 2  %Bisection Guess
                    [pout.s.v1,pout.s.v2,pout.s.xnext, pout.s.vnext] = BPswitchBisect(x,v1,v2,x1init,x2init);
            end
            if ~isempty(pout.s.v2)
                print_diag(2,'Branch Angle = %6.6f * pi\n', innerangle(pout.s.v1,pout.s.v2)/pi())
            end
            J = contjac(x);
            if ~issparse(J)
                [~,evl1]= eig(J(:,1:ndim-1)');
            else
                opt.disp = 0;
                [d,evl1] = eigs(J(:,1:ndim-1), NumEvals, 'SM', opt);  % dsb %Question here about "6" eigenvalues, why this magic number?
            end
            
            return;
        end
        
        %JH: Newton converging slowly, exit and reduce stepsize
        if ((normf > normf_old/2) || isnan(normf) ) && ~firstnewtonstep
            print_diag(3,'Locate BP: Residual Converging Too Slowly\n');
            return
        end
        
        firstnewtonstep = 0;
        %Step 9
        bord.v = vp{1}/norm(vp{1});
        bord.w = vp{2};
        
        J = contjac(x);
    end
    
    if i == MaxCorrIters
        print_diag(3,'Locate BP: Failed to find x in %d Newton iterations\n',MaxCorrIters);
        x=[];v=[];
        return;
    end
    
    i = i+1;
end


%% Branch Angle Bisection
% ---------------------------------------------------------
function [V1,V2,xnext,vnext] = BPswitchBisect(x,v1,v2,x1init,x2init)
print_diag(5,'In equilibriumL/BPswitchBisect\n');

global cds
h0 = cds.h;
%P = [v1 v2]*[v1'; v2'];

X1 = x1init;
X2 = x2init;

% calculate alpha1, alpha2 and define thetaMin = max(alpha1,alpha2)
alpha1   = innerangle(v1, [v1 v2]*([v1'; v2']*(x - X1)));
alpha2   = innerangle(v1, [v1 v2]*([v1'; v2']*(x - X2)));
thetaMin = max(alpha1,alpha2)+1e-6;
if thetaMin < (pi()/2 - 1e-6) % Tolerance
    gamma = atan(thetaMin);
else
    print_diag(3,'BPswitchBisect: Cannot locate bracketing points x1, x2 on initial curve\n')
    V1 = [];   V2 = [];
    return;
end

vright = (gamma*v2 + v1)/(sqrt(gamma^2 +1));
vleft  = (gamma*v2 - v1)/(sqrt(gamma^2 +1));
vguess = v2;

%JH:    Determines the current direction of the tangent vector
if abs(vguess(end)) < 1e-6
    Vdir = sign(sum(vguess(1:end-1)));
else
    Vdir = sign(vguess(end));
end
%JH: Ensures the tangent vector is calculated in the 'positive' direction.
if Vdir < 0
    vguess = -vguess;
end

h1 = h0;
%Tangent vector bisection calculation
while 1
    if h1 < cds.h_min/2
        print_diag(3,'BPswitchBisect: Stepsize too small\n')
        V1 = [];   V2 = [];
        return;
    end
    xguess = x + h1*vguess;
    %%%tol = 1e-4; %%% test BPswitvh #2
    %%%if abs(x(end)-xguess(end)) >= tol h1 = h1/2; continue; end;%%%
    %% x_end = x(end) %%%
    %%  xguess_end = xguess(end) %%%
    pcorr = newtcorrL(xguess, vguess, []);
    if   isempty(pcorr);   h1 = h1/2; continue; end
    xcorr = pcorr.x; vcorr = pcorr.v;
    
    theta = innerangle(v1,[v1 v2]*([v1'; v2']*(x-xcorr)));
    if theta < thetaMin && norm(vleft - vright) > 1e-6
        if v1'*[v1 v2]*([v1'; v2']*(x-xcorr)) > 0
            vright = vguess;
            print_diag(2,'Right Angle = %d\n', thetaMin - theta)%test
        else
            vleft  = vguess;
            print_diag(2,'Left  Angle = %d\n', thetaMin - theta)%test
        end
        
        vguess = vleft + vright;
        vguess = vguess/norm(vguess);
    else
        break;
    end
end
%Next point prediction.
h1 = h0;
while 1
    if h1 < cds.h_min/2
        print_diag(3,'BPswitchBisect: Cannot locate nearby point\n')
        V1 = [];   V2 = [];
        return;
    end
    xguess = x + h1*vcorr;
    pnext = newtcorrL(xguess,vcorr, []);
    if   isempty(pnext)
        h1 = h1/2;
    else
        xnext = pnext.x; vnext = pnext.v;
        theta = innerangle(v1,[v1 v2]*([v1'; v2']*(x-xnext)));
        if theta > thetaMin     %check for fallback or curves too close
            V1 = v1;  V2 = vguess; V2 = V2/norm(V2);
            return;
        else
            print_diag(3,'BPswitchBisect: Unresolved Fallback or Branches too close, Unable to switch branches.\n')
            V1 = v1;   V2 = []; xnext = []; vnext = [];
            return;
        end
        
    end
end