function hopf_for_cycle
  % N is the number of mesh points of the discretization of the system of PDE's
  N = 31; 
  % other parameters of the system of ODEs:
  L = 0.5; A = 2; B = 5.45; Dx = 0.008; Dy = 0.004;
  ode_parameters = [N L A B Dx Dy];


  % the continuation will be with respect to the second parameter:
  active_parameter = 2;

  % We specify the initial point "equilibrium".
  % This equilibrium represents a spatially homogeneous equilibrium of the PDE:
  equilibrium          = zeros(2*N,1);
  equilibrium(1:N)     = A;
  equilibrium(N+1:2*N) = B/A;
  
  odefile = @brusselator_1d;
  
  [x0,v0] = init_EP_EP_L(odefile, equilibrium, ...
                    ode_parameters, active_parameter);

  % we specify the options for the equilibrium continuation, by creating a
  % struct opts_ep_ep, which will be passed to contL.

  opts_ep_ep = contset(...
    'MaxNumPoints', 23, ...
    'Singularities', true);

  % we run the equlibrium continuation:
  contL(@equilibriumL,x0,v0, opts_ep_ep);
end