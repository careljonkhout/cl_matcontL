function bc = bratu_bc (ya, yb)
%
%% SAMPLE1_BC evaluates the boundary conditions.
%
%    Input, real YA(M), YB(M), the solution value at the left and right endpoints.
%    Output, real BC(2?), the value of the boundary conditions.
%
  bc(1) = ya(1);
  bc(2) = yb(1);
%------------------------------------------------------------------
