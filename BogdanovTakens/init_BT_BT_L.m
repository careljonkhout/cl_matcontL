function [x0, v0]= init_BT_BT_L(odefile, u, p, ap,data)
  %  Initializes a Bogdanov-Takens continuation from a BT point

  if nargin > 5
    error('init_BT_BT_L cannot be called with more than 5 arguments');
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
    end
  end
  if size(ap,2)~=3
    error(['Three active parameter are needed ' ...
           'for a Bogdanov-Takens bifurcation continuation']);
  end
  
  x0 = [u; p(ap)'];
  v0 = [];

   
  global cds
  cds = [];
  
  load_odefile(odefile); % adds fields to cds
    
  cds.P0           = p;
  cds.ActiveParams = ap;
  cds.ndim         = length(x0);
  cds.ncoo         = length(x0) - 3;
  cds.nap          = length(ap);
  cds.curve        = @bogdanovtakensL;
end