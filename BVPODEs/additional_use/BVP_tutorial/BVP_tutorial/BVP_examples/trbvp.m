function trbvp
%TRBVP  Exercise for Example 3 of the BVP tutorial.
%   This problem is studied in section 5.4 of B.A. Finlayson, The Method of
%   Weighted Residuals and Variational Principles, Academic, New York, 1972.
%   It arises when modelling a tubular reactor with axial dispersion.  An
%   isothermal situation with n-th order irreversible reaction leads to 
%   
%      y'' = Pe*(y' - R*y^n)
%   
%   Here Pe is the axial Peclet number and R is the reaction rate group.
%   The boundary conditions are
%   
%      y'(0) = Pe*(y(0) - 1),    y'(1) = 0.
%   
%   Finlayson compares results he obtains with an orthogonal collocation
%   method to results obtained by others with finite differences when 
%   Pe = 1, R = 2, and n = 2.  These results are compared here to results
%   obtained with BVP4C. BVPVAL is used to get a smoother graph of y(x).   


% Known parameter
Pe = 1;

options = bvpset('stats','on');
solinit = bvpinit(linspace(0,1,5),[0.5 0]);
sol = bvp4c(@trode,@trbc,solinit,options,Pe);

fprintf('\n');
fprintf('Other authors report y(0) = 0.63678, y(1) = 0.45759.\n');
fprintf('Values computed are  y(0) = %7.5f, y(1) = %7.5f\n',sol.y(1,1),sol.y(1,end));

clf reset
xint = linspace(0,1);
Sxint = bvpval(sol,xint);
plot(xint,Sxint(1,:))
title('Mass transfer in a tubular reactor.')
xlabel('x')
ylabel('y')
shg

% --------------------------------------------------------------------------

function dydx = trode(x,y,Pe)
%TRODE  ODE function for the exercise of Example 3 of the BVP tutorial.
dydx = [ y(2)
         Pe*(y(2) + 2*y(1)^2)];

% --------------------------------------------------------------------------

function res = trbc(ya,yb,Pe)
%TRBC  Boundary conditions for the exercise of Example 3 of the BVP tutorial.
res = [ ya(2) - Pe*(ya(1) - 1)
        yb(2) ];

