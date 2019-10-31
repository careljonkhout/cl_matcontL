% based on algorithm outlined in paragraph 6.2 of (bibtex citation follows)
% @phdthesis{lust-phd,
%	  author={Lust, Kurt},
%	  title={Numerical bifurcation analysis 
%        of periodic solutions of partial differential equations},
%	  school={K.U.Leuven},
%	  year={1997},
% }
% most variable names are derived from variable names in \cite{lust-phd}.
function x = do_one_correction(x0, x, v0)
  global cds;

  [V, reduced_jacobian, delta_q_gamma, delta_q_r, G_delta_q_r, ...
          phases_0, phi, period, active_par_val] = ...
    NewtonPicard.MultipleShooting.compute_reduced_jacobian(x);
  
  m = cds.nMeshIntervals;
  
  lhs_3_1 = zeros(1, cds.reduced_jac_size);
  col_offset = 0;
  for i = 1 : m % m == cds.nMeshIntervals
    indices_lhs_3_1          = col_offset  + ( 1 : size(V{i}, 2) );
    indices_v0               = (i - 1) * cds.nphases + (1 : cds.nphases);
    lhs_3_1(indices_lhs_3_1) = v0(indices_v0)' * V{i};
    col_offset               = col_offset + size(V{i}, 2);
  end
  
  c_n__delta_q_gamma = 0;
  
  for i = 1 : m
    indices_v0          = (i-1) * cds.nphases + (1:cds.nphases);
    c_n__delta_q_gamma  = c_n__delta_q_gamma + ...
                                  v0(indices_v0)' * delta_q_gamma(:,i);
  end
  
  left_hand_side = [
    reduced_jacobian;
    lhs_3_1    v0(end-1)  v0(end) + c_n__delta_q_gamma;
  ];


  rhs_1      = zeros(cds.reduced_jac_size, 1);
  row_offset = 0;
  for i=1:m
    ni             = next_index_in_cycle(i,m);
    indices        = row_offset + ( 1 : size(V{ni}, 2) );
    rhs_1(indices) = V{ni}' * ( phi(:,i) - phases_0(:,ni) + G_delta_q_r(:,i) );
    row_offset     = row_offset + size(V{ni}, 2);
  end
 
  
  right_hand_side = - [
    rhs_1;
    
    cds.previous_dydt_0' * ...
      (phases_0(:,1) - cds.previous_phases  + delta_q_r(:,1));
      
    v0'*(x-x0)  + v0(1:end-2)' * reshape(delta_q_r,numel(delta_q_r),1);
  ];

  delta_p_delta_T_and_delta_gamma = lsqminnorm(left_hand_side, right_hand_side);
  
  % delta_p containts the corrections in the subspace V all the mesh points on
  % the cycle. The subspace V is the subspace associated the leading eigenvalues
  % of the current approximation of the monodromy matrix of the curr. approx. of
  % the cycle. Note that we store m mesh points and have m mesh intervals.
  % Hence, after the corrections in a continuation step have converged, the
  % first mesh point does NOT coincide with the last mesh point continuation.
  % (In contrast the cl_matcontL implementation of cycle continuation by
  % orthogonal collocation, the first and the last mesh point do approximately
  % coincide, when the corrections have converged)
  delta_p     = delta_p_delta_T_and_delta_gamma(1:end-2);
  delta_T     = delta_p_delta_T_and_delta_gamma(end-1);
  delta_gamma = delta_p_delta_T_and_delta_gamma(end);

  V_delta_p = zeros(cds.nphases * m, 1);
  col_offset = 0;
  for i = 1 : m
    indices_delta_p              = col_offset + ( 1 : size(V{i}, 2) );
    indices_V_delta_p            = (i - 1) * cds.nphases + (1 : cds.nphases);
    V_delta_p(indices_V_delta_p) = V{i} * delta_p(indices_delta_p);
    col_offset                   = col_offset + size(V{i}, 2);
  end
  
  % note x = x(:) reshapes x into a column vector
  
  phases_0       = phases_0(:);
  delta_q        = delta_q_r + delta_gamma * delta_q_gamma;
  phases_0       = phases_0 + V_delta_p + delta_q(:);
  period         = period + delta_T;
  active_par_val = active_par_val + delta_gamma;
  x              = [phases_0; period; active_par_val];
end

function index = next_index_in_cycle(i,m)
  index = i + 1;
  if index == m + 1
    index = 1;
  end
end
