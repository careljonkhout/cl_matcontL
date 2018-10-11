function [x,v] = CorrectStartPoint(x0, v0)
global cds

x = [];
v = [];

% no tangent vector given, cycle through base-vectors

if cds.options.TSearchOrder==1;
i = 1;
v0 = zeros(cds.ndim,1);
while isempty(x) & i<=cds.ndim
  v0(i) = 1;
   try
  DefaultProcessor(x0,v0);
  [x,v] = newtcorr(x0, v0);
   catch
   end
  v0(i) = 0; 
  i=i+1;
end
  else
i = cds.ndim;
v0 = zeros(cds.ndim,1);
while isempty(x) & i>=1
  v0(i) = 1;
   try
  DefaultProcessor(x0,v0);
  [x,v] = newtcorr(x0, v0);
   catch
   end
  v0(i) = 0; 
  i=i-1; 
end
end