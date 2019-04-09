function v = find_tangent_vector(x)

   [~, reduced_jacobian, ~, ~, ~, ~] = ...
          PartialMonodromy.compute_reduced_jacobian(x);

%	if (rank(reduced_jacobian) < min(size(reduced_jacobian)))
%    error(['The reduced Jacobian is rank deficient. ' ...
%           'Cannot find tangent vector. '])
%  end
        
  phases__T_and_gamma = lsqminnorm(  reduced_jacobian, ...
                                     zeros( size(reduced_jacobian, 1), 1 )  );
  phases = phases__T_and_gamma(1:end-2);
  T      = phases__T_and_gamma(end-1);
  gamma  = phases__T_and_gamma(end);

  
  
  delta_phases   = phases;
  v              = [delta_phases; T; gamma];
    
