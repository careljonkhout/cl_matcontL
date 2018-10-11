function bratubvp_mark
%BRATUBVP  Exercise for Example 1 of the BVP tutorial.
%   The BVP  y'' + exp(y) = 0, y(0) = 0 = y(1) is a standard example
%   of a problem with two solutions.  It is easy enough to solve, but
%   some experimentation with the guess may be necessary to get both.


options = bvpset('stats','on');
solinit = bvpinit(linspace(0,1,5),[0.1 0]);
%%sol1 = bvp5c_mark(@bratuode,@bratubc,solinit,options); %4->5 mark
sol1 = bvp5c(@bratuode,@bratubc,solinit,options); %4->5 mark

fprintf('\n');

% Change the initial guess to converge to a different solution. 
solinit = bvpinit(linspace(0,1,5),[3 0]);
sol2 = bvp5c_mark(@bratuode,@bratubc,solinit,options); %4->5 mark

figure
plot(sol1.x,sol1.y(1,:),sol2.x,sol2.y(1,:))
title('Bratu''s equation has two solutions when \lambda = 1.')
xlabel('x')
ylabel('y')

% --------------------------------------------------------------------------

function dydx = bratuode(x,y)
%BRATUODE  ODE function for the exercise of Example 1 of the BVP tutorial.
dydx = [  y(2)
         -exp(y(1))];

% --------------------------------------------------------------------------

function res = bratubc(ya,yb)
%BRATUBC  Boundary conditions for the exercise of Example 1 of the BVP tutorial.
res = [ya(1)
       yb(1)];

