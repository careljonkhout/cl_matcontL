function [x0, v0] = init_multiple_shooting_from_hopf( odefile, x, ...
   ode_parameters, active_parameter_index, h, nMeshIntervals, subspace_size)


  global cds
  
  if length(active_parameter_index) ~= 1
      error(['One active parameter is needed for limit cycle ' ...
        'continuation using single shooting']);
  end
  

  cds.nphases = length(x) - 1;
  cds.curve = @single_shooting;
  [max_order, ~] = find_maximum_order_of_symbolic_derivatives(odefile);
  cds.options.SymDerivative = max_order;

  cds.mesh            = linspace(0, 1, nMeshIntervals + 1);
  
  handles                = feval(odefile);
  dydt_ode               = handles{2};
  jacobian_ode           = handles{3};
  
  if isnumeric(ode_parameters)
    ode_parameters = num2cell(ode_parameters);
  end
  
  A = jacobian_ode(0, x(1:end-1), ode_parameters{:});
  
  
  % calculate eigenvalues and eigenvectors
  [V,D] = eigs(A, cds.nphases);
  % find pair of complex eigenvalues
  d = diag(D);
  smallest_sum = Inf;
  for j=1:cds.nphases-1
    [val,idx] = min(abs(d(j+1:cds.nphases)+d(j)));
    if val < smallest_sum
      idx1 = j;
      idx2 = j+idx;
      smallest_sum = val;
    end
  end
  % real part? Oh dear, a neutral saddle!
  if imag(d(idx1)) == 0 && imag(d(idx2)) == 0
    x0=[];
    v0=[];
    debug('Neutral saddle\n');
    return;
  end
  % get imaginary part and corresponding eigenvector
  omega = abs(imag(d(idx2)));
  Q = V(:,idx1);

  d = real(Q)'*real(Q);
  s = imag(Q)'*imag(Q);
  r = real(Q)'*imag(Q);
  Q = Q*exp(1i*atan2(2*r,s-d)/2);
  Q = Q/norm(real(Q));

  % initial amplitude h
  % calculate initial cycle and its tangent vector
  % It is uncertain if the next 2 lines are correct, but it seems to work
  t = kron(exp(2*pi*1i*cds.mesh(1:end-1)),Q);
  
  v0 = [real(t(:));0;0];
  
  
  active_param_val = ode_parameters{active_parameter_index};
  
  x = x(1:cds.nphases);
  
  x0 = [repmat(x,length(cds.mesh)-1,1);2*pi/omega;active_param_val]+h*v0;
  
  v0 = v0/norm(v0);


  point_on_limitcycle    = x0(1:cds.nphases);
  tangent_to_limitcycle  = dydt_ode(0, point_on_limitcycle, ode_parameters{:});

  cds.nMeshIntervals  = nMeshIntervals;
  cds.probfile        = odefile;
  cds.options.PartitionMonodromy = cds.nphases > 20;
  cds.nap             = 1;
  cds.ndim            = nMeshIntervals * cds.nphases + 2;
  cds.usernorm        = [];
  cds.ncoo            = cds.nphases + 1;
  cds.ActiveParams    = active_parameter_index;
  cds.P0              = cell2mat(ode_parameters);
  cds.previous_phases = point_on_limitcycle;
  cds.previous_dydt_0 = tangent_to_limitcycle;
  cds.dydt_ode        = dydt_ode;
  cds.jacobian_ode    = jacobian_ode;
  cds.jacobian_p_ode  = handles{4};
  cds.integrator      = @ode15s;
  cds.preferred_basis_size  = subspace_size;
  cds.p               = subspace_size;
  cds.mv_count        = 0;
  cds.curve           = @multiple_shooting;
 
end

