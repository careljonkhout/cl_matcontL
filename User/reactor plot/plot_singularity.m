function plot_singularity(s, varargin)
  paramater = s.data.x(end);
  T = s.data.T;
  plot(paramater,T,'r*')
  text(paramater,T,s.label,varargin{:})
end