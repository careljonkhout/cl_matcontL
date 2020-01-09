function [x0, v0]= init_H_H_L(odefile, u, p, ap, data)
  %
  % Initializes a Hopf bifurcation continuation from a Hopf point found 
  %

  if nargin > 5
    error('init_H_H_L cannot be called with more than 5 arguments');
  end
  if isempty(data)
    error('Data from equilibriumL must be present'); 
  end

  if isempty(u)
    u = data.x(1:end-1);
  end
  if isempty(p)
    p = data.P0;
    p(data.ap) = data.x(end); % DV 2018
    if size(p, 2) > 1
      p = p';
    end
  end
  if size(ap,2)~=2
    error(['Two active parameter are needed ' ...
           'for a Hopfpoint bifurcation continuation']);
  end

  x0 = [u; p(ap); data.s.kapa];        % change continuation parameters
  v0 = [];
  
  global cds
  cds = [];

  load_odefile(odefile); % adds fields to cds

  % initialize cds
  cds.P0           = p;
  cds.ActiveParams = ap;
  cds.ndim         = length(x0);
  cds.ncoo         = length(x0) - 3;
  cds.nap          = length(ap);
  cds.curve        = @hopfL;
end