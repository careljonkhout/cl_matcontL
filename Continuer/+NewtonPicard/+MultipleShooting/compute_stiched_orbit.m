% computes the orbit of the cycle by "stiching together" all the partial orbits
% started from all the mesh-points.

function compute_stiched_orbit(x)
  global cds contopts;
   int_opt = odeset( ...
    'AbsTol', contopts.integration_abs_tol, ...
    'RelTol', contopts.integration_rel_tol  ...
  );
  [phases_0, period, parameters] = ...
    NewtonPicard.MultipleShooting.extract_phases_period_and_parameters(x);
  dydt       = @(t,x) cds.dydt_ode(t, x, parameters{:});
  t_cycle = [];
  y_cycle = [];
  indices = 1:cds.nphases;
  for i=1:cds.nMeshPoints
    time_interval = period * [cds.mesh(i) cds.mesh(i+1)];
    [t_mesh_interval, y_mesh_interval] = ...
      cds.integrator(dydt, time_interval, phases_0(indices), int_opt);
    if i < cds.nMeshPoints
      t_cycle = [t_cycle; t_mesh_interval(1:end-1,:)]; %#ok<AGROW>
      y_cycle = [y_cycle; y_mesh_interval(1:end-1,:)]; %#ok<AGROW>
    else
      t_cycle = [t_cycle; t_mesh_interval(1:end  ,:)]; %#ok<AGROW>
      y_cycle = [y_cycle; y_mesh_interval(1:end  ,:)]; %#ok<AGROW>
    end
    % We do not know the end size of t_cycle and y_cycle beforehand, therefore,
    % we ignore the 'variable appears to change size on every loop iteration' -
    % warning.
    indices = indices + cds.nphases;
  end
  
  cds.t_cycle = t_cycle;
  cds.y_cycle = y_cycle;
end