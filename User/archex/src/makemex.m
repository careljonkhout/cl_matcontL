function makemex
if ~exist('femex')
  mex femex.c fem1.c gaussquad.c
end


