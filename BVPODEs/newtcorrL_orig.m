function [x,v,i,A,funcnorm] = newtcorrL(x0, v0)
%
% Newton corrections, internal routine
%
global cds contopts

MaxNewtonIters = contopts.Cont_MaxNewtonIters;
MaxCorrIters   = contopts.Cont_MaxCorrIters;
Solver         = contopts.Cont_Solver;
VarTolerance   = contopts.Cont_VarTolerance;
FunTolerance   = contopts.Cont_FunTolerance;
sparse_flag = contopts.CIS_SparseSolvers;

x = x0;
v = v0;

R = [];
R(cds.ndim) = 1;

Q = [feval(cds.curve_func, x); 0];      % orig %MP 7-2017
%%% Q = [feval(cds.curve_func, x); 0];    % orig %MP 4-2016
%%Q1 = feval(cds.curve_func, x);          %MP 4-2016
%%   if ~sparse_flag                      %MP 4-2016
%%   Q = [Q1 0];                          %MP 4-2016
%%   Q=Q';                                %MP 4-2016  
%%elseif sparse_flag                     %MP 4-2016
%%   Q = [Q1; 0];                         %MP 4-2016
%%end                                     %MP 4-2016

funcnorm = unorm(Q,2);

% WM: for-loops are faster then while-loops
for i = 1:MaxCorrIters
    
    if i <= MaxNewtonIters
        A = contjac(x);
        if isempty(A)
            B = [v'];
        else
            B = [A; v(1:cds.ndim)'];
        end
    end
    
    % repeat twice with same jacobian, calculating
    % the jacobian is usually a lot more expensive
    % than solving a system
    for t=1:2
        if isempty(A) || length(Q) == 1
            x = [];
            v = [];
            return;
        end
        
        if strcmp(Solver,'Moore-Penrose')
            
            lastwarn('');
            D = bordCIS1(B(1:cds.ndim,1:cds.ndim),[Q(1:cds.ndim) R'],1);        % dx = B\Q; dv = B\R' CIS;
            %%D = B\[Q R'];                     % test   ????
            if ~isempty(lastwarn)
                x = [];
                v = [];
                return;
            end
            
            v = D(:,2);
            v = v/norm(v);
            
            dx = D(:,1);
        else
            dx = B\Q;
        end
        
%         % DSB -- line search to ensure *some* decrease (if not Armijo)
%         s = 1;
%         decrease = 0;
%         for jj = 1:10
%             xtrial = x - s*dx;
%             Q = [feval(cds.curve_func, xtrial); 0];
%             if length(Q) > 1
%                 decrease = 1;
%                 funcnormt = unorm(Q,2);
%                 if funcnormt < funcnorm
%                     break; 
%                 end
%             end
%             s = s/2;
%         end       
%         if decrease
%             funcnorm = funcnormt;
%             x = xtrial;
%         else
%             x = [];
%             v = [];
%             return;
%         end
        
        
        
        
        %if norm(dx) < cds.options.VarTolerance & ...       %MF
        %      norm(Q) < cds.options.FunTolerance           %MF
        %if normU(dx,1) < cds.options.VarTolerance & ...    %MF %JH
        %   normU(Q(1:end-1),2) < cds.options.FunTolerance  %MF %JH
        
        % *** User Norm check ***                               %JH
        varnorm = unorm(dx,1);
        % *** end user norm check ***
        
        %%if varnorm < VarTolerance && funcnorm < FunTolerance;
        if varnorm < VarTolerance & funcnorm < FunTolerance;
            A = contjac(x);
            if isempty(A)
                x = [];
                v = [];
                return;
            end
            
            v = bordCIS1([A;v'],R',1);
            v = v/norm(v);
            return;
        end
        
        x = x - dx;

Q = [feval(cds.curve_func, x); 0];      % orig %MP 7-2017
%%% Q = [feval(cds.curve_func, x); 0];    % orig %MP 4-2016
%%Q1 = feval(cds.curve_func, x);          %MP 4-2016
%%   if ~sparse_flag                      %MP 4-2016
%%   Q = [Q1 0];                          %MP 4-2016
%%   Q=Q';                                %MP 4-2016  
%%elseif sparse_flag                     %MP 4-2016
%%   Q = [Q1; 0];                         %MP 4-2016
%%end                                     %MP 4-2016

        if length(Q) == 1 % curve_func failed
            x = []; v = [];
            return
        end
        
        funcnorm = unorm(Q,2);
        
    end
end

funcnorm
varnorm
x = [];
v = [];
A = [];


function normv = unorm(v, which)
global cds
if ~isempty(cds.usernorm)
    normv = normU(v,which);
else
    normv = norm(v);
end
