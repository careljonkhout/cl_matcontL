function point = do_corrections(x0,v)
  
  global cds contopts
 
  curve_function_norm = max(abs(feval(cds.curve_func,x0)));
  x = x0;
  corrections = 0;
  if curve_function_norm > 10
    point = [];
    return;
  end
  while (~ (curve_function_norm < contopts.FunTolerance)) ...
      && corrections < contopts.MaxCorrIters
    period = x(end-1);
    fprintf('function_norm: %.8e period: %.8e corrections: %d\n', ...
      curve_function_norm, period, corrections);
    corrections = corrections + 1;
    if     isequal(cds.curve, @single_shooting)
      x = NewtonPicard.SingleShooting.do_one_correction(x0,x,v);
    elseif isequal(cds.curve, @multiple_shooting)
      x = NewtonPicard.MultipleShooting.do_one_correction(x0,x,v);
    else
      print_diag(0,'Newton_Picard_Corrections: wrong curvefile.\n');
      point = [];
      return
    end
    if isempty(x)
      break;
    end

    if period < 0
      fprintf('period less than zero, correction aborted\n');
      break
    end
    new_curve_function_norm = max(abs(feval(cds.curve_func,x)));
    if (new_curve_function_norm > 2 * curve_function_norm) 
      print_diag(0,'function_norm: %.8e period: %.8e corrections: %d\n', ...
        new_curve_function_norm, period, corrections);
      print_diag(0,['Curve function norm is strongly increasing.' ...
        ' Aborting Corrections.\n']);
      break;
    end
    curve_function_norm = new_curve_function_norm;
  end
  if curve_function_norm < contopts.FunTolerance
    point.R = curve_function_norm; 
    point.x = x;
    point.v = v;
    point.iters = corrections;
    if isequal(cds.curve, @single_shooting)


      active_par_val               = x(end);
      period                       = x(end-1);
      phases_0                     = x(1:end-2);
      parameters                   = cds.P0;
      parameters(cds.ActiveParams) = active_par_val;
      parameters                   = num2cell(parameters);

      integration_opt = odeset(...
        'AbsTol',      1e-10,    ...
        'RelTol',      1e-10,    ...
        'BDF',         'off',   ...
        'MaxOrder',     5,      ...
        'NormControl',  'off',  ...
        'Refine',       1,      ...
        'Jacobian',     @(t,y) feval(cds.jacobian_ode,t,y,parameters{:}) ...
      );


      cds.cycle_orbit = cds.integrator(...
        @(t, y) cds.dydt_ode(t, y, parameters{:}), ...
        linspace(0, period, cds.nDiscretizationPoints), ...
        phases_0, integration_opt);

      %compute_subspace(period, parameters);
      %disp(cds.eigenvalues)
    end
  else
    point = [];
    return
  end
end