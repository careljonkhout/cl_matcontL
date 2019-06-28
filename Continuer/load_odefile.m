function load_odefile(odefile)

  global cds
  handles = feval(odefile);

  dydt            = handles{2};
  jacobian        = handles{3};
  jacobian_params = handles{4};
  hessians        = handles{5};
  hessians_params = handles{6};
  der3            = handles{7};
  der4            = handles{8};
  der5            = handles{9};
  user_norm       = handles{10};
  user_func = cell(length(handles)-10,1);
  for ii = 11:length(handles)
      user_func{ii-10} = handles{ii};
  end

  % initialize cds
  cds.probfile     = odefile;
  cds.func         = dydt;
  cds.Jacobian     = jacobian;
  cds.JacobianP    = jacobian_params;
  cds.Hessian      = hessians;
  cds.HessianP     = hessians_params;
  cds.der3         = der3;
  cds.der4         = der4;
  cds.der5         = der5;
  cds.usernorm     = user_norm;
  cds.userf        = user_func;
  cds.curve        = @equilibriumL;
end


