odefile = @nonadiabatic_tubular_reactor;

load('/home/carel/Documents/cl_matcontL/User/na_tubular_reactor/Data/na_tubular_reactor_orb_lc_04-Apr-2019_20_03_53.mat')
%init_LPC_LPC_L(odefile, x, s, ap, ntst, ncol,varargin)
  [x,v] = init_LPC_LPC_L(odefile, s.data.x, s, [1 2], 40, 4);
  contL(@limitpointcycle,x,v)