% based on algorithm outlined in paragraph 6.2 of (bibtex citation follows)
% @phdthesis{lust-phd,
%	  author={Lust, Kurt},
%	  title={Numerical bifurcation analysis 
%        of periodic solutions of partial differential equations},
%	  school={K.U.Leuven},
%	  year={1997},
% }
% most variable names are derived from variable names in \cite{lust-phd}.
function v = find_tangent_vector(x)
  global cds
  m = cds.nMeshIntervals;
  
  [V, reduced_jacobian, delta_q_gamma, delta_q_r, ~, ~, ~, ~, ~] = ...
    NewtonPicard.MultipleShooting.compute_reduced_jacobian(x);
  
  basis_size = size(V,2);

  delta_p__delta_T_and_delta_gamma = null(full(reduced_jacobian));
  delta_p     = delta_p__delta_T_and_delta_gamma(1:end-2);
  delta_T     = delta_p__delta_T_and_delta_gamma(end-1);
  delta_gamma = delta_p__delta_T_and_delta_gamma(end);

  V_delta_p = zeros(cds.nphases*m,1);
  for i=1:m % m == cds.nMeshIntervals
    indices1 = (i-1) * cds.nphases + (1:cds.nphases);
    indices2 = (i-1) * basis_size  + (1:basis_size );
    V_delta_p(indices1) = V(:,:,i) * delta_p(indices2);
  end
  
  delta_q         = delta_q_r + delta_gamma * delta_q_gamma;
  delta_phases        = V_delta_p + reshape(delta_q,numel(delta_q),1);
  v = [delta_phases; delta_T; delta_gamma];
end
    
