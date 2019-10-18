% inputs (Lust)
% based on algorithm 6.2 on pages 197 and 198 of (bibtex citation follows)
% @phdthesis{lust-phd,
%	  author={Lust, Kurt},
%	  title={Numerical bifurcation analysis
%          of periodic solutions of partial differential equations},
%	  school={K.U.Leuven},
%	  year={1997},
% }
%
% inputs:
% - V: bases for the subspace of the monodromy matrix at point x_i on the
% cycle. That is, V(:,:,i) contains column vectors that span the subspace.
%
% - rhs: the right hand side of the Q system
% note: the minus in rhs is implicit, i.e. this function solves .... = - rhs
%
% - delta_t: the lengths of the time mesh intervals. 
%
% - parameters: cell array of the values of the parameters of the system of ODEs

function [delta_q, G_delta_q] = solve_Q_system(V, rhs, delta_t, parameters)
  global cds contopts;
  
  m = cds.nMeshIntervals;
  for i=1:m
    ni = next_index_in_cycle(i,m);
    rhs(:,i) = rhs(:,i) - V(:,:,ni) * V(:,:,ni)' * rhs(:,i);
  end
  
  G_delta_q  = zeros(cds.nphases,m);
  delta_q    = zeros(cds.nphases,m);
  
  minimum_residual = Inf;
  for iteration_number = 1:contopts.MaxPicardIterations
   
    for i=2:m %  m == cds.nMeshIntervals;
      delta_q(:,i) = G_delta_q(:,i-1) + rhs(:,i-1);
      delta_q(:,i) = delta_q(:,i) - V(:,:,i) * V(:,:,i)' * delta_q(:,i);
      G_delta_q(:,i) = NewtonPicard.MultipleShooting.monodromy_map(i, ...
                                      delta_q(:,i), delta_t(i), parameters);
    end
    condensed_residual = G_delta_q(:,m) + rhs(:,m);
    condensed_residual = condensed_residual ...
                         - V(:,:,1) * V(:,:,1)' * condensed_residual;
    residual_norm = max(max(abs(delta_q(:,1) - condensed_residual)));
    if residual_norm < contopts.PicardTolerance
      break
    elseif residual_norm < minimum_residual
      minimum_residual = residual_norm;
      delta_q(:,1) = condensed_residual;
      G_delta_q(:,1) = NewtonPicard.MultipleShooting.monodromy_map( ...
        1, delta_q(:,1), delta_t(1), parameters);
    else % if residual_norm >= minumim_residual
      print_diag(4,'poor convergence in solve_q_system. residual: %.5e\n', ...
                    residual_norm);
	    break
    end
  end
  print_diag(5,'did %d iterations in solve_q_system\n',iteration_number);
end
  
  
function index = next_index_in_cycle(i,m)
  index = i + 1;
  if index == m + 1
    index = 1;
  end
end
    
