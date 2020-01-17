% based on page 283 of (bibtex citation follows)
% @incollection{lust2000,
%	title={Computation and bifurcation analysis of periodic solutions of large-scale systems},
%	author={Lust, Kurt and Roose, Dirk},
% booktitle={Numerical methods for bifurcation problems and large-scale dynamical systems},
%	pages={265--301},
%	year={2000},
%	publisher={Springer}
% }
% this function is currently not used (october 2019), since recomputing the
% subspaces at each point seems to be faster.
function V = NP_SS_continue_subspace(period, parameters)
  global cds contopts
  mv_count = 0;
  V_extended = [cds.V ...
    rand(cds.n_phases, cds.p + 2 - size(cds.V,2))];
  %V_extended = [rand(cds.n_phases, cds.p + cds.p_extra)];
  p = cds.p;
  
  int_opt = odeset(...
    'AbsTol',       contopts.integration_abs_tol,    ...
    'RelTol',       contopts.integration_rel_tol,    ...
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
     
      V_extended(:,p_eff+1:end) = orth(W(:,p_eff+1:end));
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
        if (strcmp(error.identifier,'MATLAB:remoteparfor:AllParforWorkersAborted'))
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
    
    if ~ contopts.contL_ParallelComputing || parfor_failed
      for i=p_eff+1:size(V_extended,2)
        [~,orbit] = cds.integrator(...
          dydt_monodromy_map, [0 period], V_extended(:,i), int_opt);
        W(:,i) = orbit(end,:)';
        mv_count = mv_count +1;
      end     
    end
        
      
    U = V_extended'*W;
    [Y,S] = schur(U);
    % Order schur vectors find the permutation whose inverse sorts the Schur
    % vectors by absolute value of the associated eigenvalues in descending
    % order:
    [~, permutation] = sort(abs(ordeig(S)));
    % Invert permutation:
    permutation(permutation) = 1:length(permutation);
    % Apply the permutation to sort the Schur vector in the required order:
    [Y,S] = ordschur(Y,S,permutation);
    V_extended = V_extended * Y;
    W = W * Y;
    k = size(V_extended,2) - 1;
    % note: svds(A,1) is a built-in matlab function that computes the largest
    % singular value of A
    while k >= 1 ...
        && (S(k+1,k) ~= 0 || S(k+1,k) == 0 ...
        && svds(W(:,1:k) - V_extended(:,1:k) * S(1:k,1:k), 1) ...
              >= contopts.NewtonPicardBasisTolerance )
      k = k - 1;

    end
    fprintf('subspace iteration %d mv_count: %d\n',iteration,mv_count);
    p_eff = k;
    %fprintf('p_eff: %d subspace norm at k=p_eff: %.10f\n', p_eff, svds(W(:,1:p_eff) - V_extended(:,1:p_eff) * S(1:p_eff,1:p_eff), 1));
  end
  V = V_extended(:,1:p_eff);
  
  
