function plot_singularity(s, varargin)
  L = s.data.parametervalues(3);
  T = s.data.T;
  plot(L,T,'r*')
  text(L,T,s.label,varargin{:})
end