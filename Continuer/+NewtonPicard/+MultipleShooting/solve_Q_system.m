  % inputs (Lust)
% based on algorithm 6.2 on pages 197 and 198 of (bibtex citation follows)
% @phdthesis{lust-phd,
%	  author={Lust, Kurt},
%	  title={Numerical bifurcation analysis
%          of periodic solutions of partial differential equations},
%	  school={K.U.Leuven},
%	  year={1997},
% }
% inputs:
% - starting values delta_q_i, i=0..m-1, each in \R^N (N is number of spatial 
%   dimensions of ode
% - optionally G_i(delta_q_i) i=0..m-1 each in \R^n
% - convergence thresholds eps_i 
% - bases V for projectors
% - right hand sides rhs
% - routine "monodromy_map" to compute G delta_q_i
%
% P will be the projectors onto the small subspaces
%

function [delta_q, G_delta_q] = ...
             solve_Q_system(V, rhs, partial_period, parameters)
  global cds contopts;
  m = cds.nMeshPoints;
  for i=1:m
    ni = next_index_in_cycle(i,m);
    rhs(:,i) = rhs(:,i) - V(:,:,ni) * V(:,:,ni)' * rhs(:,i);
  end
  
%   if max(max(abs(rhs))) < eps
%     delta_q   = zeros(cds.nphases,m);
%     G_delta_q = zeros(cds.nphases,m);
%     return
%   end
  G_delta_q  = zeros(cds.nphases,m);
  delta_q    = zeros(cds.nphases,m);

  for iteration_number = 1:15
   
    for i=2:m
      delta_q(:,i) = G_delta_q(:,i-1) + rhs(:,i-1);
      delta_q(:,i) = delta_q(:,i) - V(:,:,i) * V(:,:,i)' * delta_q(:,i);
      G_delta_q(:,i) = NewtonPicard.MultipleShooting.monodromy_map(i, ...
        delta_q(:,i), partial_period, parameters);
    end
    condensed_residual = G_delta_q(:,m) + rhs(:,m);
    condensed_residual = condensed_residual ...
                         - V(:,:,1) * V(:,:,1)' * condensed_residual;
    if max(max(abs(delta_q(:,1) - condensed_residual))) ...
        < contopts.PicardTolerance
      break
    else
      delta_q(:,1) = condensed_residual;
      G_delta_q(:,1) = NewtonPicard.MultipleShooting.monodromy_map( ...
        1, delta_q(:,1), partial_period, parameters);
    end
  end
  print_diag(5,'did %d iterations in solve_q_systems\n',iteration_number);
end
  
  
function index = next_index_in_cycle(i,m)
  index = i + 1;
  if index == m + 1
    index = 1;
  end
end

function index = previous_index_in_cycle(i,m) %#ok<DEFNU>
  index = i - 1;
  if index == 0
    index = m;
  end
end
  
    
