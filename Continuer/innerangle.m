function theta = innerangle(x,y)

xunit = x/norm(x);
yunit = y/norm(y);

theta = acos(xunit'*yunit);
if theta > pi/2
    theta = abs(theta - pi);
end