h=makeGrid(100,10,3,1e-5);
x = cumsum(h)


function h = makeGrid(N, L, n, hmin)
    x = linspace(0, 1, N+1);
    x = (L - hmin*N)*x.^n;
    h = hmin + diff(x);
end

function dpdx = dY(h1, h2, Y0, Y1, Y2)
    a2 = (h1 + 2*h2)/(3*h2*(h1 + h2));
    a1 = (h2 - h1)/(3*h1*h2);
    a0 = -(2*h1 + h2)/(3*h1*(h1 + h2));
    dpdx = a2 * Y2 + a1 * Y1 + a0 * Y0;
end

function ddpdx = ddY(h1, h2, Y0, Y1, Y2)
    a2 = 2/(h2 * (h1 + h2));
    a1 = -2/(h1 * h2);
    a0 = 2/(h1 * (h1 + h2));
    ddpdx = a2 * Y2 + a1 * Y1 + a0 * Y0;
end

function dpdx = leftdY(h1, h2, Y0, Y1, Y2)
    a2 = -h1 / (h2 * (h1 + h2));
    a1 = (h1 + h2) / (h1 * h2);
    a0 = -(2*h1 + h2) / (h1 * (h1 + h2));
    dpdx = a2 * Y2 + a1 * Y1 + a0 * Y0;
end

function Y0 = leftdYs(h1, h2, Y1, Y2, dY0)
    a2 = -h1 / (h2 * (h1 + h2));
    a1 = (h1 + h2) / (h1 * h2);
    a0 = -(2*h1 + h2) / (h1 * (h1 + h2));
    Y0 = (dY0 - a2 * Y2 - a1 * Y1) / a0;
end

function Y0 = leftdYl(h1, h2, Y1, Y2, lambda)
    a2 = -h1 / (h2 * (h1 + h2));
    a1 = (h1 + h2) / (h1 * h2);
    a0 = -(2*h1 + h2) / (h1 * (h1 + h2));
    Y0 = -(a2 * Y2 + a1 * Y1) / (a0 - 1 / lambda);
end

function Y0 = leftddYz(h1, h2, Y1, Y2)
    a2 = 2/(h2 * (h1 + h2));
    a1 = -2/(h1 * h2);
    a0 = 2/(h1 * (h1 + h2));
    Y0 = -(a2 * Y2 + a1 * Y1) / a0;
end

function dpdx = rightdY(h1, h2, Y0, Y1, Y2)
    a2 = (2*h2 + h1) / (h2 * (h1 + h2));
    a1 = -(h1 + h2) / (h1 * h2);
    a0 = h2 / (h1 * (h1 + h2));
    dpdx = a2 * Y2 + a1 * Y1 + a0 * Y0;
end

function YN = rightdYs(h1, h2, Y0, Y1, dYN)
    a2 = (2*h2 + h1) / (h2 * (h1 + h2));
    a1 = -(h1 + h2) / (h1 * h2);
    a0 = h2 / (h1 * (h1 + h2));
    YN = (dYN - a1 * Y1 - a0 * Y0) / a2;
end

function YN = rightdYsl(h1, h2, Y0, Y1, dYN, lambda)
    a2 = (2*h2 + h1) / (h2 * (h1 + h2));
    a1 = -(h1 + h2) / (h1 * h2);
    a0 = h2 / (h1 * (h1 + h2));
    YN = (dYN - a1 * Y1 - a0 * Y0) / (a2 - 1/lambda);
end