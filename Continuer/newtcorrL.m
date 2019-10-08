% Newton corrections, internal routine
% if the Newton corrections fail to converge, an empty array will be returned
function pout = newtcorrL(x0, v0, CISdata)

global cds contopts

pout = [];
MaxNewtonIters = contopts.MaxNewtonIters;
MaxCorrIters   = contopts.MaxCorrIters;
MoorePenrose   = contopts.MoorePenrose;
VarTolerance   = contopts.VarTolerance;
FunTolerance   = contopts.FunTolerance;

x = x0;
v = v0;

R = [];
R(cds.ndim) = 1;

if cds.newtcorrL_needs_CISdata
    % if CIS data needed for evaluation, update CIS data as well
    CISdata = feval(cds.curve_CIS_step, x, CISdata);
    if isempty(CISdata)
        print_diag(3, 'newtcorrL: CIS step failed')
        return;
    end
end

% if ~contopts.CIS_UsingCIS                    % MP  2018
% Q = [feval(cds.curve_func, x); 0];
% else
Q = [feval(cds.curve_func, x, CISdata); 0];
% end                                          % MP 2018

funcnorm = normU(Q);
smallest_funcnorm_so_far = funcnorm;

for i = 1:MaxCorrIters
    if i > 1
      format_string = ['newton iteration %d function norm 10^%1.2f ' ...
                                         'correction norm 10^%1.2f\n'];
      print_diag(5,format_string, i, log10(funcnorm), log10(varnorm))
    else
      print_diag(5,'newton iteration %d function norm 10^%1.2f\n', ...
                                                            i, log10(funcnorm))
    end
    if i <= MaxNewtonIters
        A = contjac(x, CISdata);
        if isempty(A)
            B = v';
        else
           % B = [A; v(1:cds.ndim)'];  MP why ??????????
            B = [A; v'];              % MP 2018
        end
    end
    
    % repeat twice with same jacobian, calculating
    % the jacobian is usually a lot more expensive
    % than solving a system

    for t=1:2
        if isempty(A) || length(Q) == 1
            pout = [];
            return;
        end
        
        if MoorePenrose

            lastwarn('');
            D = bordCIS1(B(1:cds.ndim,1:cds.ndim),[Q(1:cds.ndim) R'],1);       
            % dx = B\Q; dv = B\R' CIS;
         
            if ~isempty(lastwarn)
                pout = [];
                return;
            end
            
            v = D(:,2);
            v = v/norm(v);
            
            dx = D(:,1);
        else
          if isequal(cds.curve, @limitcycleL)
            dx = linear_solver_collocation(B,Q);
          else
            dx = B \ Q;
          end
        end
        
        %         % DSB -- line search to ensure *some* decrease (if not Armijo)
        %         s = 1;
        %         decrease = 0;
        %         for jj = 1:10
        %             xtrial = x - s*dx;
        %             Q = [feval(cds.curve_func, xtrial); 0];
        %             if length(Q) > 1
        %                 decrease = 1;
        %                 funcnormt = normU(Q,2);
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
        varnorm = normU(dx);
        % *** end user norm check ***
        
    
        if varnorm < VarTolerance && funcnorm < FunTolerance
            A = contjac(x, CISdata);
            if isempty(A)
                pout = [];
                return;
            end
            if isequal(cds.curve,@limitcycleL)
              v = linear_solver_collocation([A;v'],R');
            else
              v = bordCIS1([A;v'],R',1);
            end
            if contopts.newtcorrL_use_max_norm
              v = v/max(abs(v));
            else
              v = v/norm(v);
            end
            pout.x = x;
%             newtcorrL_1 =1
%            length_x =  length(x)
            pout.v = v;
            pout.iters = i;
            pout.R = funcnorm;
            pout.tvals = [];
            pout.uvals = [];
            return;
        end
        
        x = x - dx;
        if cds.newtcorrL_needs_CISdata % DV: if CIS data needed for evaluation, update CIS data as well
            CISdata = feval(cds.curve_CIS_step, x, CISdata);
            if isempty(CISdata)
                print_diag(3, 'newtcorrL: CIS step failed')
                return;
            end
        end
        
        Q = [feval(cds.curve_func, x, CISdata); 0];
        
        if length(Q) == 1 % curve_func failed
            pout = [];
            return
        end
        
        funcnorm = normU(Q);
        smallest_funcnorm_so_far = min(funcnorm, smallest_funcnorm_so_far);
        if funcnorm > contopts.max_rel_funcnorm_increase ...
                               * smallest_funcnorm_so_far
	        print_diag(1,[ ...
              'Current curve function norm is now %d times ' ...
              'the lowest curve function norm that ' ...
              'was attained in this continuation step. ' ...
              'Aborting Corrections.\n'], contopts.max_rel_funcnorm_increase);
          
          format_string = 'function norm 10^%1.2f correction norm 10^%1.2f\n';
          print_diag(1,format_string, log10(funcnorm), log10(varnorm))
          pout = [];
          return;
        end
        
        if isnan(funcnorm)
          print_diag(1,'function norm is NaN. Aborting corrections\n');
          pout = [];
          return;
        end
        
    end
end