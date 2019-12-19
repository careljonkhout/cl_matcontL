function multipliers = NP_SS_compute_multipliers(x, nMults_to_compute)
  global cds contopts;
  print_diag(3,'computing multipliers\n');
  [phases, period, parameters] = NP_SS_extract_phases_period_and_parameters(x);
  
  integration_opt = odeset(...
    'AbsTol',      contopts.multipliers_abs_tol,    ...
    'RelTol',      contopts.multipliers_rel_tol     ...
  );
  if ~ isempty(cds.jacobian_ode)
    integration_opt = odeset(integration_opt, ...
          'Jacobian',     @(t,y) feval(cds.jacobian_ode,t,y,parameters{:}));
  end
  
  if ~ cds.using_cvode
    cds.cycle_orbit = cds.integrator(...
         @(t, y) cds.dydt_ode(t, y, parameters{:}), ...
         [0 period], ...
         phases(:,1), ...
         integration_opt);
  end
 
  cds.phases_0 = x(1:end-2);
  
  monodromy_map = @(x) NP_SS_monodromy_map( ...
                        x, period, parameters, ...
                        contopts.multipliers_abs_tol, ...
                        contopts.multipliers_rel_tol);
                        
  
  nMults_to_compute = min(nMults_to_compute, length(cds.phases_0)/2);
  
   [~, multiplier_matrix] = eigs(monodromy_map, cds.n_phases, ...
                                                nMults_to_compute);
                                              %'largestabs', ... 
                                              %'Tolerance', 1e-12);
  multipliers = diag(multiplier_matrix);
  
  
  multipliers = sort(multipliers, 'descend');
  multipliers = multipliers(1:nMults_to_compute);
  
  
  
  %options.k = n_eigenvalues;
  %options.blsz = 1;
%  [Q, T] = ahbschur(monodromy_map, cds.n_phases, options);
 
%   while j <= nMults_to_compute
%     multipliers(j) = T(j, j)
%     j = j + 1;
%     if T(j, j-1) ~= 0
%       a = T(j-1,j-1);
%       b = T(j-1, j);
%       c = T(j, j-1);
%       d = T(j, j);
%       multipliers(j-1) =  1/2 (-sqrt(a^2 - 2*a*d + 4*b*c + d^2) + a + d)
%       multipliers(j) =    1/2 (+sqrt(a^2 - 2*a*d + 4*b*c + d^2) + a + d)
%       j = j + 1
%     end
%   end
  
  multipliers = multipliers(1:nMults_to_compute);
  
  print_diag(1, multipliers2str(multipliers));
end