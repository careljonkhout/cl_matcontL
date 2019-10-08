function x = init_single_shooting_branching(odefile, bpc, h, ...
                                                        time_integration_method)
  parameters            = bpc.P0;
  parameters(bpc.ap)    = bpc.x(end);
  parameters            = num2cell(parameters);
  point_on_limitcycle   = bpc.x(1:end-2);
  ode_handles           = feval(odefile);
  dydt_ode              = ode_handles{2};
  tangent_to_limitcycle = feval(dydt_ode,0, point_on_limitcycle, parameters{:});
  using_cvode           = endsWith(func2str(time_integration_method), 'cvode');
  nphases               = length(bpc.x) - 2;
  basis_size            = max(nphases, min(40, nphases/2));
  
  global cds;
  cds = [];
  cds.nphases                    = length(bpc.x) - 2;
  cds.probfile                   = odefile;
  cds.options.PartitionMonodromy = cds.nphases > 20;
  cds.nap                        = 1;
  cds.ndim                       = cds.nphases + 2;
  cds.usernorm                   = [];
  cds.ncoo                       = cds.nphases + 1;
  cds.ActiveParams               = bpc.ap;
  cds.P0                         = cell2mat(parameters);
  cds.previous_phases            = point_on_limitcycle(:);
  cds.previous_dydt_0            = tangent_to_limitcycle(:);
  cds.dydt_ode                   = dydt_ode;
  cds.jacobian_ode               = ode_handles{3};
  cds.jacobian_p_ode             = ode_handles{4};
  cds.integrator                 = time_integration_method;
  cds.preferred_basis_size       = basis_size;
  cds.p                          = basis_size;
  cds.mv_count                   = 0;
  cds.curve                      = @single_shooting;
  cds.using_cvode                = using_cvode;
  
  
  v = NewtonPicard.SingleShooting.find_branching_vector(bpc.x, bpc.v);
  v = v / max(v);
  x = bpc.x + h * v;
  
  
end