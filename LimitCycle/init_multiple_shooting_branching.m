function x = init_multiple_shooting_branching(odefile, bpc, h, ...
                                                        time_integration_method)
  parameters            = bpc.P0;
  parameters(bpc.ap)    = bpc.x(end);
  parameters            = num2cell(parameters);
  n_phases              = (length(bpc.x) - 2) / bpc.n_mesh_intervals;
  point_on_limitcycle   = bpc.x(1:n_phases);
  ode_handles           = feval(odefile);
  dydt_ode              = ode_handles{2};
  tangent_to_limitcycle = bpc.v(1:n_phases);
  using_cvode           = endsWith(func2str(time_integration_method), 'cvode');
  basis_size            = min(2 * length(bpc.multipliers), n_nphases);
  
  global cds;
  cds = [];
  cds.n_phases                    = n_phases;
  cds.probfile                   = odefile;
  cds.options.PartitionMonodromy = n_phases > 20;
  cds.nap                        = 1;
  cds.ndim                       = n_phases * bpc.n_mesh_intervals + 2;
  cds.usernorm                   = [];
  cds.ncoo                       = n_phases * bpc.n_mesh_intervals + 1;
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
  cds.curve                      = @multiple_shooting;
  cds.using_cvode                = using_cvode;
  cds.n_mesh_intervals           = bpc.n_mesh_intervals;
  cds.new_n_mesh_intervals       = bpc.n_mesh_intervals;
  cds.mesh                       = bpc.time_mesh;
  
  
  v = NP_MS_find_branching_vector(bpc.x, bpc.v);
  v = v / max(v);
  x = bpc.x + h * v;
  
  
end