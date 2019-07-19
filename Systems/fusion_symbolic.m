function out = fusion()

out{1} = @init;
out{2} = @nonuniformGrid;
out{3} = [];
out{4} = [];
out{5} = [];
out{6} = [];
out{7} = [];
out{8} = [];
out{9} = [];
out{10} = [];%@usernorm;
out{11} = @user1;
out{12} = @user2;
out{13} = @user3;

end

% --------------------------------------------------------------------------
function x0 = init( N, Gamma_inf, q_inf, D0, D1, D2, a, b, zeta1, mu1, epsilon, ZS, gamma1, lambdan, lambdaT, cn, cT)

L = 10;
Zguess = 1;

h = makeGrid(N, L, 3, 1e-5);
x = cumsum(h);

Tinf = (gamma1 - 1)* q_inf / Gamma_inf;
T0 = Tinf / (1 + lambdan / (lambdaT * zeta1));
cg = (zeta1*cT - cn) / (1 + zeta1 * lambdaT/lambdan);

D   = @(z) D0 + D1*tanh(z) + D2 * (tanh(z)).^2;
G   = @(z) a - b*(z - ZS) - (z - ZS).^3;
f =@(z) -Tinf / (Gamma_inf * lambdan^2) * (cn + cg) .* D(z) + G(z);
Z0 = fzero(f, Zguess);

De = D(Z0);
Ge = G(Z0);

% Coefficients of Taylor expansion of 1/D(Z)
Dt0 = D(ZS);
Dt1 = 4*D1*(exp(ZS))^2/((exp(ZS))^2+1)^2+8*D2*(exp(ZS))^2*((exp(ZS))^2-1)/((exp(ZS))^2+1)^3;
Dt2 = -4*D1*(exp(ZS))^2*((exp(ZS))^2-1)/((exp(ZS))^2+1)^3+D2*(-8*((exp(ZS))^2-1)^2*(exp(ZS))^2/((exp(ZS))^2+1)^4+16*(exp(ZS))^4/((exp(ZS))^2+1)^4);
Dt3 = (8/3)*D1*(exp(ZS))^2*((exp(ZS))^4-4*(exp(ZS))^2+1)/((exp(ZS))^2+1)^4+D2*((16/3)*((exp(ZS))^2-1)*(exp(ZS))^2*((exp(ZS))^4-4*(exp(ZS))^2+1)/((exp(ZS))^2+1)^5-32*(exp(ZS))^4*((exp(ZS))^2-1)/((exp(ZS))^2+1)^5);

Dt = [Dt3, Dt2, Dt1, Dt0];
DT = De * Ge * Dt;
GT = [-1 0 -b a];
p = GT - DT;

r = roots(p);
r = sort(r);

if ~isreal(r)
    ZL = fzero(G, -1);
    Z = ZL * ones(1, length(x));
else

ZL = r(1);
ind = find(abs(r - Z0) == min(abs(r - Z0)));
if ind == 1
    Z = ZL * ones(1, length(x));
else
    if ind == 2
        Z2 = r(2);
        Z3 = r(3);
    elseif ind == 3
        Z2 = r(3);
        Z3 = r(2);
    end
    z2 = Z2 - ZL;
    z3 = Z3 - ZL;
    c = sqrt(-p(1)/ (1 * mu1));

    k = c * sqrt(z2 * z3);
    qp = sqrt(2*z3 - z2) * ( sqrt(2) * sqrt(2*z3 - z2) + sqrt(3*z3) ) / ((sqrt(2*z3) + sqrt(z2)) * (sqrt(z3) + sqrt(2*z2)));
    qm = sqrt(2*z3 - z2) * ( sqrt(2) * sqrt(2*z3 - z2) - sqrt(3*z3) ) / ((sqrt(2*z3) + sqrt(z2)) * (sqrt(z3) + sqrt(2*z2)));
    Z = r(1) + (r(3) - r(1)) * (qm + 1) ./ (qm + exp(k*x)) * (qp + 1) ./ (qp + exp(-k*x));
end
end

n0 = -Gamma_inf * lambdan / D(Z0);
n = n0 - Gamma_inf * cumsum(1./D(Z) .* h);
T = Tinf + (T0 - Tinf) * (n / n0).^(-zeta1);
U = n .* T / (gamma1 - 1);

x0(1:3:3*(N-1), 1) = n(1:end-1); 
x0(2:3:3*(N-1), 1) = U(1:end-1);
x0(3:3:3*(N-1), 1) = Z(1:end-1);

f =@(t,y) nonuniformGrid(t, y , N, Gamma_inf, q_inf, D0, D1, D2, a, b, zeta1, mu1, epsilon, ZS, gamma1, lambdan, lambdaT, cn, cT);
[~, Y] = ode15s(f, [0 100], x0);
x0 = Y(end, :)';
end

% --------------------------------------------------------------------------
function [ dydt ] = nonuniformGrid(~, y , N, Gamma_inf, q_inf, D0, D1, D2, a, b, zeta, mu, epsilon, ZS, gamma, lambdan, lambdaT, cn, cT)
% unused parameter is t

L = 10;
h = makeGrid(N, L, 3, 1e-5);
M = N-1;

D   = @(z) D0 + D1*tanh(z) + D2*(tanh(z)).^2;
chi = @(z) D(z) / (zeta * (gamma - 1));
G   = @(z) a - b*(z - ZS) - (z - ZS).^3;

dDdx  = @(z) D1*(1 - tanh(z)^2) + 2*D2*tanh(z)*(1 - tanh(z)^2);
dchidx= @(z) dDdx(z) / (zeta * (gamma - 1));

% split
n = y(1:3:3*M);
U = y(2:3:3*M);
T = (gamma - 1) * U ./ n;
Z = y(3:3:3*M);

dndt =  sym(zeros(N-1, 1));
dUdt =  sym(zeros(N-1, 1));
dZdt =  sym(zeros(N-1, 1));

% boundary conditions x = 0
h1 = h(1);
h2 = h(2);
Zleft = leftddYz(h1, h2, Z(1), Z(2));
nleft =  leftdYl(h1, h2, n(1), n(2), lambdan);
Tleft =  leftdYl(h1, h2, T(1), T(2), lambdaT);
%Uleft = nleft * Tleft / (gamma - 1);

% boundary conditions x = L
Zright =  rightdYs(h(M), h(M+1), Z(M-1), Z(M), 0);
nright =  rightdYs(h(M), h(M+1), n(M-1), n(M), -Gamma_inf / D(Zright));
Tright = rightdYsl(h(M), h(M+1), T(M-1), T(M), -q_inf / (chi(Zright) * nright), ...
    (gamma - 1) * chi(Zright) * nright / Gamma_inf );

% extend with boundary points
n = [nleft  n  nright];
T = [Tleft  T  Tright];
Z = [Zleft  Z  Zright];

% internal points
for i = 2:(M+1)
    h1 = h(i-1);
    h2 = h(i);
    
    dN  =  dY(h1, h2, n(i-1), n(i), n(i+1));
    ddN = ddY(h1, h2, n(i-1), n(i), n(i+1));
    
    dT  =  dY(h1, h2, T(i-1), T(i), T(i+1));
    ddT = ddY(h1, h2, T(i-1), T(i), T(i+1));
    
    dZ  =  dY(h1, h2, Z(i-1), Z(i), Z(i+1));
    ddZ = ddY(h1, h2, Z(i-1), Z(i), Z(i+1));
    
    dD  =  dDdx(Z(i)) * dZ;
    dchi=  dchidx(Z(i)) * dZ;
    
    dndt(i-1) = dN * dD + D(Z(i)) * ddN;
    dUdt(i-1) = chi(Z(i)) * n(i) * ddT + (dchi * n(i) + chi(Z(i)) * dN) * dT + ...
                (1 / (gamma - 1)) * ( D(Z(i)) * T(i) * ddN + (dD * T(i) + D(Z(i)) * dT) * dN );
    dZdt(i-1) = mu * ddZ + cn * T(i) / n(i)^2 * dN + cT / n(i) * dT + G(Z(i));
end

dydt(1:3:3*M) = dndt;
dydt(2:3:3*M) = dUdt;
dydt(3:3:3*M) = dZdt / epsilon;

end
% --------------------------------------------------------------------------

% -------------------------------------------------------------------------
function normuser = usernorm(arg)
n = length(arg);
normuser = sqrt(1/n)*norm(arg);
end
% -------------------------------------------------------------------------
function uf1 = user1(t,x,N, L, Gamma_inf, q_inf, D0, D1, a, b, zeta, mu, epsilon, ZS, gamma, lambdan, lambdaT, lambdaZ, cn, cT)
uf1 = q_inf - 0.7;
end

function uf1 = user2(t,x,N, L, Gamma_inf, q_inf, D0, D1, a, b, zeta, mu, epsilon, ZS, gamma, lambdan, lambdaT, lambdaZ, cn, cT)
uf1 = q_inf - 0.8;
end

function uf1 = user3(t,x,N, L, Gamma_inf, q_inf, D0, D1, a, b, zeta, mu, epsilon, ZS, gamma, lambdan, lambdaT, lambdaZ, cn, cT)
uf1 = q_inf - 0.9;
end
% ---------------------------------------------------------------------
% ----

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