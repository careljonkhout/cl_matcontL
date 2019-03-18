% based on algorithm outlined in paragraph 6.2 of (bibtex citation follows)
% @phdthesis{lust-phd,
%	  author={Lust, Kurt},
%	  title={Numerical bifurcation analysis 
%        of periodic solutions of partial differential equations},
%	  school={K.U.Leuven},
%	  year={1997},
% }
% most variable names are derived from variable names in \cite{lust-phd}.
function x = do_one_correction(x0,x,v0)

  [V, reduced_jacobian, delta_q_gamma, delta_q_r, G_delta_q_r] = ...
    compute_reduced_jacobian(x);
  lhs_3_1 = zeros(1,basis_size * m);
  for i=1:m % m == cds.nShootingPoints
    indices_lhs_3_1 = (i-1)*basis_size  + (1:basis_size);
    indices_v0      = (i-1)*cds.nphases + (1:cds.nphases);
    lhs_3_1(indices_lhs_3_1) = v0(indices_v0)' * V(:,:,i);
  end
  
  left_hand_side = [
    reduced_jacobian;
    lhs_3_1    v0(end-1)         v0(end)     ;
  ];


  rhs_1 = zeros(basis_size*m,1);
  for i=1:m
    indices        = (i-1) * basis_size + (1:basis_size);
    ni             = next_index_in_cycle(i,m);
    rhs_1(indices) = V(:,:,ni)'*(phi(:,i) - phases_0(:,ni) + G_delta_q_r(:,i));
  end
 
  
  right_hand_side = - [
    rhs_1;
    
    cds.previous_dydt_0' * ...
      (phases_0(:,1) - cds.previous_phases  + delta_q_r(:,1));
      
    v0'*(x-x0)  + v0(1:end-2)' * reshape(delta_q_r,numel(delta_q_r),1);
  ];

  delta_p__delta_T_and_delta_gamma = left_hand_side \ right_hand_side;
  delta_p     = delta_p__delta_T_and_delta_gamma(1:end-2);
  delta_T     = delta_p__delta_T_and_delta_gamma(end-1);
  delta_gamma = delta_p__delta_T_and_delta_gamma(end);

  V_delta_p = zeros(cds.nphases*m,1);
  for i=1:m
    indices1 = (i-1) * cds.nphases + (1:cds.nphases);
    indices2 = (i-1) * basis_size  + (1:basis_size );
    V_delta_p(indices1) = V(:,:,i) * delta_p(indices2);
  end
  phases_0 = reshape(phases_0,numel(phases_0),1);
  delta_q         = delta_q_r + delta_gamma * delta_q_gamma;
  phases_0        = phases_0 + V_delta_p + reshape(delta_q,numel(delta_q),1);
  period          = period + delta_T;
  active_par_val  = active_par_val + delta_gamma;
  x = [phases_0; period; active_par_val];
  %v = find_tangent_vector(phases_0, period, parameters, V);
end

% extracts 
% - y ( the current approximation of points on the cycle )
% - period
% - parameters ( of the ode system in which cycles are continued )
% from the continuation state vector x.
% The non-active parameters, i.e. the parameters that do not change during
% the continuation are extracted from the global struct cds
% (i.e.) curve description structure.
% The parameters are returned as a cell array, so that that can be passed to
% cds.dydt_ode in an syntactically elegant manner.

function [y,period,parameters] = getComponents(x)
  global cds
  y                            = x(1:cds.nphases*cds.nShootingPoints);
  y                            = reshape(y,cds.nphases,cds.nShootingPoints);
  period                       = x(end-1);
  parameter_value              = x(end);
  parameters                   = cds.P0;
  parameters(cds.ActiveParams) = parameter_value;
  parameters                   = num2cell(parameters);
end


    
function x_end = shoot(x, period, parameters)
  global cds
  f =@(t, y) cds.dydt_ode(t, y, parameters{:});
  integration_opt = odeset(...
    'AbsTol',      1e-10,    ...
    'RelTol',      1e-10,    ...
    'BDF',         'off',   ...
    'MaxOrder',     5,      ...
    'NormControl',  'off',  ...
    'Refine',       1,      ...
    'Jacobian',     @(t,y) feval(cds.jacobian_ode,t,y,parameters{:}) ...
  );
  [~, trajectory] = ode15s(f, [0 period], x, integration_opt);
  x_end = trajectory(end,:)';
end

 
function dphidp = d_phi_d_gamma(x, partial_period, parameters)
  global cds
  ap = cds.ActiveParams;
  h = 1e-6;
  parameters{ap} = parameters{ap} - h;
  phi_1 = shoot(x, partial_period, parameters);
  parameters{ap} = parameters{ap} + 2*h;
  phi_2 = shoot(x, partial_period, parameters);
  dphidp = (phi_2 - phi_1)/h/2;
end

  
function V = compute_subspace(i, period, parameters)
  global cds
  
  [eigenvectors, eigenvalues, no_convergence] = eigs( ...
    @(x) NewtonPicard.MultipleShooting.monodromy_map( ...
       i, x, period, parameters), ...
    cds.nphases, ...
    cds.p + 1);


  if no_convergence
    V = [];
    fprintf(['Newton_Picard_Correction.m:', ...
      ' eigenvalues of monodromy matrix did not converge.\n'])
    return
  end

  eigenvalues = diag(eigenvalues);
  basis = zeros(cds.nphases, cds.p + 1);

  i = 1;

  while i <= cds.p
    basis(:, i) = real(eigenvectors(:,i));
    i = i + 1;
    if max(abs(imag(eigenvalues(i-1)))) > 1e-14
      basis(:, i) = imag(eigenvectors(:,i-1));
      i = i + 1;
    end
  end
  V = orth(basis(:,1:i-1));
end


function index = next_index_in_cycle(i,m)
  index = i + 1;
  if index == m + 1
    index = 1;
  end
end

function index = previous_index_in_cycle(i,m) %#ok<DEFNU>
  index = i - 1;
  if index == 0
    index = m;
  end
end
  
    
