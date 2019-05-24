function theta = innerangle(x,y)

global cds

if isequal(cds.curve, @limitcycleL)
  x = x(end-1,end);
  y = y(end-1,end);
end

xunit = x/norm(x);
yunit = y/norm(y);

theta = acos(xunit'*yunit);
if theta > pi/2
    theta = abs(theta - pi);
end