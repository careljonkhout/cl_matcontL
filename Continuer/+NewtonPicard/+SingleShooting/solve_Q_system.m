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
  %delta_q    = zeros(size(V,1),1);
  M_delta_q  = zeros(size(V,1),1);
  %residual   = rhs + M_delta_q - delta_q;
  %residual   = residual - V * V' *residual;
  %if max(abs(residual)) < contopts.PicardTolerance
  %  return;
  %end
  for iteration_number = 1:contopts.MaxPicardIterations
    delta_q = M_delta_q + rhs;
    delta_q = delta_q - V*V'*delta_q;
    M_delta_q = NewtonPicard.SingleShooting.monodromy_map( ...
        delta_q, period, parameters);
    residual = rhs + M_delta_q - delta_q;
    residual = residual - V*V'*residual;
    print_diag(5,'q_system residual %.5e\n', max(abs(residual)));
    if max(abs(residual)) < contopts.PicardTolerance
      break
    end
  end
  print_diag(5,'did %d iterations in solve_q_systems\n',iteration_number);
end