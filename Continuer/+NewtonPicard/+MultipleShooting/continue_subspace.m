% based on page 283 of (bibtex citation follows)
% @incollection{lust2000,
%	title={Computation and bifurcation analysis of periodic solutions of large-scale systems},
%	author={Lust, Kurt and Roose, Dirk},
% booktitle={Numerical methods for bifurcation problems and large-scale dynamical systems},
%	pages={265--301},
%	year={2000},
%	publisher={Springer}
% }
function V = continue_subspace(i,phases,period, parameters)
  global cds contopts
  V_extended = [cds.V(:,:,i) ...
    rand(cds.nphases, cds.p + 2 - size(cds.V(:,:,i),2))];
  %V_extended = [rand(cds.nphases, cds.p + cds.p_extra)];
  p = cds.p;
  
  integration_opt = odeset(...
    'AbsTol',      1e-10,    ...
    'RelTol',      1e-10,    ...
    'BDF',         'off',   ...
    'MaxOrder',     5,      ...
    'NormControl',  'off',  ...
    'Refine',       1,      ...
    'Jacobian',     @(t,y) feval(cds.jacobian_ode,t,y,parameters{:}) ...
  );
  

  
  cycle_trajectory = cds.integrator(...
      @(t, y) cds.dydt_ode(t, y, parameters{:}), ...
      linspace(0, period, cds.nDiscretizationPoints), ...
      phases, integration_opt);
    
  
  int_opt = odeset(...
    'AbsTol',       1e-10,    ...
    'RelTol',       1e-10,    ...
    'BDF',          'off',   ...
    'MaxOrder',     5,      ...
    'NormControl',  'off',  ...
    'Refine',       1, ...
    'Jacobian',     @(t,y) feval(cds.jacobian_ode, ...
                      t, deval(cycle_trajectory,t), parameters{:}) ...
  );                
  dydt_monodromy_map = @(t, y) ...
    cds.jacobian_ode(t, deval(cycle_trajectory,t), parameters{:}) * y;

  
  p_eff = 0;
  W = V_extended;
  iteration = 0;
  while p_eff < p && iteration < contopts.NewtonPicardMaxSubspaceIterations
    iteration = iteration + 1;
    if  iteration > 3
      fprintf('subspace iteration %d\n',iteration);
      V_extended(:,p_eff+1:end) = orth(W(:,p_eff+1:end));
      %V_extended = orth(W);
    end
    
    try 
      for i=p_eff+1:size(V_extended,2)
        % The function monodromy_map cannot be used here, since it depends on
        % the global variable cds, and global variables are not copied so the
        % the workspace of the workers that parfor uses.
        [~,trajectory] = cds.integrator(...
          dydt_monodromy_map, [0 period], V_extended(:,i), int_opt);
        W(:,i) = trajectory(end,:)';
      end
    catch error
      if (strcmp(error.identifier,'MATLAB:remoteparfor:AllParforWorkersAborted'))
        % Something went wrong with the parfor workers.
        % We try again with ordinary for.
        fprintf('Parfor aborted, retrying with ordinary for.\n');
        for i=p_eff+1:size(V_extended,2)
          [~,trajectory] = cds.integrator(...
            dydt_monodromy_map, [0 period], V_extended(:,i), int_opt);
          W(:,i) = trajectory(end,:)';
        end
      else
        % in case of some other error, we want to know about it
        rethrow(error)
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
    % note: svds(A,1) is a built-in matlab function
    % that computes the largest singular value of A    
    while k >= 1 ...
        && (S(k+1,k) ~= 0 || S(k+1,k) == 0 ...
        && svds(W(:,1:k) - V_extended(:,1:k) * S(1:k,1:k), 1) >= contopts.NewtonPicardBasisTolerance )
      k = k - 1;

    end
    p_eff = k;
    %fprintf('p_eff: %d subspace norm at k=p_eff: %.10f\n', p_eff, svds(W(:,1:p_eff) - V_extended(:,1:p_eff) * S(1:p_eff,1:p_eff), 1));
  end
  V = V_extended(:,1:p_eff);
  
  
