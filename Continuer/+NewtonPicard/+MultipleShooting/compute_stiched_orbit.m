% computes the orbit of the cycle by "stiching together" all the partial orbits
% started from all the mesh-points.

function compute_stiched_orbit(x, abs_tol, rel_tol)
  global cds contopts;
  
  if nargin == 3
    abs_tol = contopts.integration_abs_tol;
    rel_tol = contopts.integration_rel_tol;
  elseif nargin ~= 1
    error( ['The number of input arguments to ' ...
      'NewtonPicard.MultipleShooting.compute_stiched_orbit is not correct.\n' ...
      'The number of input arguments should be either 1 or 3.\n' ...
      'The actual number of input arguments is %d.\n'], nargin);
  end
  
  int_opt = odeset( ...
    'AbsTol', abs_tol, ...
    'RelTol', rel_tol  ...
  );
  [phases_0, period, parameters] = ...
    NewtonPicard.MultipleShooting.extract_phases_period_and_parameters(x);
  dydt       = @(t,x) cds.dydt_ode(t, x, parameters{:});
  t_cycle = [];
  y_cycle = [];
  indices = 1:cds.nphases;
  for i=1:cds.nMeshIntervals
    time_interval = period * [cds.mesh(i) cds.mesh(i+1)];
    [t_mesh_interval, y_mesh_interval] = ...
      cds.integrator(dydt, time_interval, phases_0(indices), int_opt);
    if i < cds.nMeshIntervals
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