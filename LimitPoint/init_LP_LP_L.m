function [x0,v0]= init_LP_LP_L(odefile, u, p, ap, data)
  %  Initializes a Limit Point / Fold continuation from a LP point

  if nargin > 5        
      error('init_LP_LP_L cannot be called with more than 5 arguments');
  end
  if isempty(u)
      if nargin < 5 % no data struct available
          error('No initial point found');
      else
          u = data.x(1:end-length(data.ap));
      end
  end
  if isempty(p)
      if nargin < 5 % no data struct available
          error('No parameters found');
      else
          p = data.P0;
          p(data.ap) = data.x(end); % MP 2020
      end       
  end                                                        
  if size(ap,2)~=2
      error(['Two active parameter are needed ' ...
             'for a Limitpoint bifurcation continuation']);
  end
  x0 = [u; p(ap)'];        % change continuation parameters
  v0 = [];
  
  global cds
  cds = [];

  load_odefile(odefile); % adds fields to cds

  % initialize cds
  cds.P0           = p;
  cds.ActiveParams = ap;
  cds.ndim         = length(x0);
  cds.ncoo         = length(x0) - 2;
  cds.nap          = length(ap);
  cds.curve        = @limitpointL;

  % cds.lastwh       = data.lastwh;  DV: We do not need this anymore, because
  % we can use the CIS algorithm to assure smooth borders
  cds.BranchParams = [];   % DV: branch switching on LP curve not yet tested