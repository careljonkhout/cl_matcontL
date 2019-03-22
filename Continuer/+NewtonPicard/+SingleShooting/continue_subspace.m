% based on page 283 of (bibtex citation follows)
% @incollection{lust2000,
%	title={Computation and bifurcation analysis of periodic solutions of large-scale systems},
%	author={Lust, Kurt and Roose, Dirk},
% booktitle={Numerical methods for bifurcation problems and large-scale dynamical systems},
%	pages={265--301},
%	year={2000},
%	publisher={Springer}
% }
function V = continue_subspace(period, parameters)
  global cds contopts
  V_extended = [cds.V rand(cds.nphases, cds.p + cds.p_extra - size(cds.V,2))];
  %V_extended = [rand(cds.nphases, cds.p + cds.p_extra)];
  p = cds.p;
  int_opt = odeset(...
    'AbsTol',       contopts.continue_subspace_abs_tol,    ...
    'RelTol',       contopts.continue_subspace_rel_tol,    ...
    'Jacobian',     @(t,y) feval(cds.jacobian_ode, ...
                      t, deval(cds.cycle_orbit,t), parameters{:}) ...
  );                
  dydt_monodromy_map = @(t, y) ...
    cds.jacobian_ode(t, deval(cds.cycle_orbit,t), parameters{:}) * y;

  
  p_eff = 0;
  W = V_extended;
  iteration = 0;
  while p_eff < p && iteration < contopts.NewtonPicardMaxSubspaceIterations
    iteration = iteration + 1;
    if  iteration > 3
      fprintf('subspace iteration %d\n',iteration);
      W_unlocked = orth(W(:,p_eff+1:end));
      V_extended(:,p_eff+1:p_eff+1+size(W_unlocked,2)-1) = W_unlocked;
      %V_extended = orth(W);
    end
    integrator = cds.integrator;
    if contopts.contL_ParallelComputing
      try 
        parfor i=p_eff+1:size(V_extended,2)
          % The function monodromy_map cannot be used here, since it depends on
          % the global variable cds, and global variables are not copied so the
          % the workspace of the workers that parfor uses.
          [~,orbit] = feval(integrator, ...
            dydt_monodromy_map, [0 period], V_extended(:,i), int_opt);
          W(:,i) = orbit(end,:)';

        end
        parfor_failed = false;
      catch error
        if (strcmp(error.identifier, ...
            'MATLAB:remoteparfor:AllParforWorkersAborted'))
          % Something went wrong with the parfor workers.
          % We try again with ordinary for.
          fprintf('Parfor aborted, retrying with ordinary for.\n');
          parfor_failed = true;
        else
          % in case of some other error, we want to know about it
          rethrow(error)
        end
      end
    end
    if (~ contopts.contL_ParallelComputing) || parfor_failed
      for i=p_eff+1:size(V_extended,2)
        [~,orbit] = integrator(dydt_monodromy_map, ...
          [0 period], V_extended(:,i), int_opt);
        
        W(:,i) = orbit(end,:)';
      end
    end
        
      
    U = V_extended'*W;
    [Y,S] = schur(U);
    % order schur vectors
    [~,I] = sort(abs(ordeig(S)));
    I(I) = 1:length(I);
    [Y,S] = ordschur(Y,S,I);
    V_extended = V_extended * Y;
    W = W * Y;
    k = size(V_extended,2) - 1;
 
    while k ~= 0
      % note: svds(A,1) is a built-in matlab function
      % that computes the largest singular value of A   
      basis_norm = svds(W(:,1:k) - V_extended(:,1:k) * S(1:k,1:k), 1);
      print_diag(3,'k: %d subspace norm at k: %.10f\n', ...
        k, basis_norm);
      if basis_norm < contopts.NewtonPicardBasisTolerance 
        break
      end
      k = k - 1;
    end
    p_eff = k;
    print_diag(3,'p_eff: %d subspace norm at k=p_eff: %.10f\n', ...
      p_eff, basis_norm);
  end
  V_extended_orth = orth(V_extended);
  V = V_extended_orth(:,1:p_eff);
  
  
