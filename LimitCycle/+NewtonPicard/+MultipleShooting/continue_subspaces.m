% based on page 208 of Kurt Lust PhD thesis
% this function is currently not used (october 2019), since recomputing the
% subspaces at each point seems to be faster.
function V = continue_subspaces(delta_t, parameters)
  global cds contopts
  m = cds.n_mesh_intervals;         
  extended_basis_size = cds.p + 4;
  V_extended = zeros([cds.n_phases extended_basis_size m+1]);
  V_extended(:,:,1) = [cds.V(:,:,1) ...
    rand(cds.n_phases, extended_basis_size - size(cds.V,2))];
  
  p_eff = 0;
  iteration = 0;
  
  while p_eff < cds.p && iteration < contopts.NewtonPicardMaxSubspaceIterations
    iteration = iteration + 1;
    if  iteration > 3
      print_diag(2,'subspace iteration %d\n',iteration);
    end
    
    V_extended_1_orth = orth(V_extended(:,:,1));
    for i = 1:100
      if size(V_extended_1_orth,2) ~= extended_basis_size
        V_extended_1_orth = [V_extended_1_orth ...
          rand(size(V_extended_1_orth,1), ...
          cds.p + 2 - size(V_extended_1_orth,2))]; %#ok<AGROW>
        % This loop will rarely run at all. It will only happen if the rank of
        % the basis drops by accident. Therefore we ignore the 'array size will
        % grow on each iteration warning'.        
        V_extended_1_orth = orth(V_extended_1_orth);
      else
        break
      end
      if i==10
        error(['Internal error in file %s.\n.' ...
          '  Current iteration should not run 100 times'], ...
          mfilename('fullpath'));
      end
    end
    
    V_extended(:,:,1) = V_extended_1_orth;

    integrator = cds.integrator;
    
    for i=1:m % m == cds.n_mesh_intervals
      int_opt = odeset(...
        'AbsTol',       contopts.integration_abs_tol,    ...
        'RelTol',       contopts.integration_rel_tol,    ...
        'Jacobian',     @(t,y) feval(cds.jacobian_ode, ...
                          t, deval(cds.orbits(i),t), parameters{:}) ...
      );
  
      dydt_monodromy_map = @(t, y) ...
        cds.jacobian_ode(t, deval(cds.orbits(i),t), parameters{:}) * y;
            
      for j = p_eff+1:size(V_extended, 2)
        V_extended(:,j , i+1) = monodromy_map(integrator, ...
          dydt_monodromy_map, V_extended(:,j,i), delta_t(i), int_opt);
      end
    end
      
    U = V_extended(:,:,1)' * V_extended(:,:,m+1);
    [Y,S] = schur(U);
    % order schur vectors
    [~,I] = sort(abs(ordeig(S)));
    I(I) = 1:length(I);
    [Y,S] = ordschur(Y,S,I);
    V_extended(:,:,1  ) = V_extended(:,:,1  ) * Y;
    V_extended(:,:,m+1) = V_extended(:,:,m+1) * Y;
    k = size(V_extended,2) - 1;
    % note: svds(A,1) is a built-in matlab function
    % that computes the largest singular value of A    
    while k ~= 0
      basis_norm = svds( ...
        V_extended(:,1:k,m+1) - V_extended(:,1:k,1) * S(1:k,1:k), 1);
      print_diag(3,'k: %d basis_norm: %.10f\n', k, basis_norm);
      if basis_norm < contopts.NewtonPicardBasisTolerance
        break
      end
      k = k - 1;
    end
    p_eff = k;
    print_diag(3,'cssps: p_eff: %d subspace norm at k=p_eff: %.10f\n', ...
      p_eff, basis_norm);
    V_extended(:,:,1) = V_extended(:,:,m+1);
  end
  V = zeros(cds.n_phases,p_eff,m);
  for i=1:m % m == cds.n_mesh_intervals
    V_i_orth = orth(V_extended(:,:,i));
    V(:,:,i) = V_i_orth(:,1:p_eff);
  end
end
  
function My = monodromy_map(integrator,dydt,y, time_interval, int_opt)
  [~, monodromy_orbit] = feval(integrator, dydt, [0 time_interval], y, int_opt);
  My = monodromy_orbit(end,:)';
end
  
  
