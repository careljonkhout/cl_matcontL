% based on algorithm outlined in paragraph 6.2 of (bibtex citation follows)
% @phdthesis{lust-phd,
%	  author={Lust, Kurt},
%	  title={Numerical bifurcation analysis 
%        of periodic solutions of partial differential equations},
%	  school={K.U.Leuven},
%	  year={1997},
% }
% most variable names are derived from variable names in \cite{lust-phd}.
function x = do_one_correction(x0,x,v0)
  global cds;
  
  [V, reduced_jacobian, delta_q_gamma, delta_q_r, G_delta_q_r, ...
          phases_0, phi, period, active_par_val] = ...
    NewtonPicard.MultipleShooting.compute_reduced_jacobian(x);
  basis_size = size(V,2);
  m = cds.nShootingPoints;
  
  lhs_3_1 = zeros(1,basis_size * m);
  for i=1:m % m == cds.nShootingPoints
    indices_lhs_3_1 = (i-1)*basis_size  + (1:basis_size);
    indices_v0      = (i-1)*cds.nphases + (1:cds.nphases);
    lhs_3_1(indices_lhs_3_1) = v0(indices_v0)' * V(:,:,i);
  end
  
  left_hand_side = [
    reduced_jacobian;
    lhs_3_1    v0(end-1)         v0(end)     ;
  ];


  rhs_1 = zeros(basis_size*m,1);
  for i=1:m
    indices        = (i-1) * basis_size + (1:basis_size);
    ni             = next_index_in_cycle(i,m);
    rhs_1(indices) = V(:,:,ni)'*(phi(:,i) - phases_0(:,ni) + G_delta_q_r(:,i));
  end
 
  
  right_hand_side = - [
    rhs_1;
    
    cds.previous_dydt_0' * ...
      (phases_0(:,1) - cds.previous_phases  + delta_q_r(:,1));
      
    v0'*(x-x0)  + v0(1:end-2)' * reshape(delta_q_r,numel(delta_q_r),1);
  ];

  delta_p__delta_T_and_delta_gamma = left_hand_side \ right_hand_side;
  delta_p     = delta_p__delta_T_and_delta_gamma(1:end-2);
  delta_T     = delta_p__delta_T_and_delta_gamma(end-1);
  delta_gamma = delta_p__delta_T_and_delta_gamma(end);

  V_delta_p = zeros(cds.nphases*m,1);
  for i=1:m
    indices1 = (i-1) * cds.nphases + (1:cds.nphases);
    indices2 = (i-1) * basis_size  + (1:basis_size );
    V_delta_p(indices1) = V(:,:,i) * delta_p(indices2);
  end
  phases_0 = reshape(phases_0,numel(phases_0),1);
  delta_q         = delta_q_r + delta_gamma * delta_q_gamma;
  phases_0        = phases_0 + V_delta_p + reshape(delta_q,numel(delta_q),1);
  period          = period + delta_T;
  active_par_val  = active_par_val + delta_gamma;
  x = [phases_0; period; active_par_val];
  %v = find_tangent_vector(phases_0, period, parameters, V);
end

function index = next_index_in_cycle(i,m)
  index = i + 1;
  if index == m + 1
    index = 1;
  end
end
