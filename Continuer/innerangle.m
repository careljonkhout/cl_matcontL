function theta = innerangle(x,y)

global cds



if isequal(cds.curve, @limitcycleL)
  x = [x(end-1) x(end)]';
  y = [y(end-1) y(end)]';
end

% convert x and y to column vectors if they are not column vectors already

x = x(:);
y = y(:);

xunit = x/norm(x);
yunit = y/norm(y);

theta = acos(xunit'*yunit);
if theta > pi/2
    theta = abs(theta - pi);
end