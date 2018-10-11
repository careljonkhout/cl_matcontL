
function z = bordCIS1(B,r,m)
% 
% W. Govaerts, Numerical methods for bifurcations of 
% dynamical equilibria, Siam 2000, p. 63, Algorithm BED
% solve B*z = r with block 2-by-2 matrix B, (1,1) block nonsingular

f = r(1:end-m,:);
g = r(end-m+1:end,:);
A = B(1:end-m,1:end-m);
b = B(1:end-m,end-m+1:end);
c = B(end-m+1:end,1:end-m)';
d = B(end-m+1:end,end-m+1:end);
P      = colamd(A);
[L, U] = lu(A(:,P));            % A(:,P) = L*U, A(:,P)' = U'*L' 

w = L'\(U'\c(P,:));             % Solve A(:,P)'*w = c(P,:)
deltas = d - w'*b;
y = deltas\(g - w'*f);
x(P,:) = U\(L\(f - b*y));       % Solve A(:,P)*x(P,:) = f - b*y

z = [x; y];
