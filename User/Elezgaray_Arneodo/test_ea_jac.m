handles = feval(@elezgaray_arneodo);
dydt = handles{2};
jac = handles{3};
global N;
N = 20;

y = [ -2 * ones(N,1) ; -4 * ones(N,1)];

D = 0.14;
eps = 0.01;
alpha = 0.01;

p = {D, eps, alpha};


dy = feval(dydt,0,y,D,eps,alpha);
my_jac = feval(jac,0,y,D,eps,alpha);

% spy(my_jac)

h = 1e-6;

nphase = size(y,1);
for i = 1:nphase
    x1 = y; x1(i) = x1(i)-h;
    x2 = y; x2(i) = x2(i)+h;
    jac_fin_diff(:,i) = feval(dydt, 0, x2, p{:})-feval(dydt, 0, x1,  p{:});
end
jac_fin_diff = jac_fin_diff/(2*h);

spy(jac_fin_diff)