% inputs (Lust)
% based on algorithm 3.1 on pages 106 and 107 of (bibtex citation follows)
% @phdthesis{lust-phd,
%	  author={Lust, Kurt},
%	  title={Numerical bifurcation analysis
%          of periodic solutions of partial differential equations},
%	  school={K.U.Leuven},
%	  year={1997},
% }
%
function [delta_q, M_delta_q] = ...
             solve_Q_system(V, rhs, period, parameters)
  global contopts;
  
  residual   = rhs - V * V' *rhs;
  if max(abs(residual)) < contopts.PicardTolerance
    print_diag(4,'did no iterations in solve_q_system\n');
    return;
  end
  
  
  M_delta_q  = zeros(size(V,1),1);
  minimal_residual = Inf;
  for iteration_number = 1:contopts.MaxPicardIterations
    delta_q = M_delta_q + rhs;
    delta_q = delta_q - V*V'*delta_q;
    M_delta_q = NewtonPicard.SingleShooting.monodromy_map( ...
        delta_q, period, parameters);
    residual = rhs + M_delta_q - delta_q;
    residual = residual - V*V'*residual;
    residual_norm = max(abs(residual));
    print_diag(5,'q_system residual %.5e\n', residual_norm);
    if residual_norm < contopts.PicardTolerance
      break
    elseif residual_norm < minimal_residual
      minimal_residual = residual_norm;
    else
      print_diag(4,'poor convergence in solve_q_system. residual: %.5e\n', ...
                    residual_norm);
	    break
    end
  end
  print_diag(4,'did %d iterations in solve_q_system\n',iteration_number);
end