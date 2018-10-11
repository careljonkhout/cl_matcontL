function [x,fail] = bordCIS2(B,r,m)
% 
% W. Govaerts & Price, 1993, Algorithm BEMW
% solve B*x = r with block 2-by-2 matrix B, (1,1) block singular

fail = 0;

i = 1;
A{1}   = B(1:end-m,1:end-m);
P      = colamd(A{1});
[L, U] = lu(A{1}(:,P));     % A(:,P) = L*U, A(:,P)' = U'*L

M{i}  = B(1:end-m+i,1:end-m+i);
b{i}  = M{i}(1:end-1,end);
bs{i} = M{i}(end,1:end-1)';
d{i}  = M{i}(end,end)';
f     = r;

vs{i} = L'\(U'\bs{i}(P,:));      % Solve A(:,P)'*vs = bs(P,:), Step 1
deltas{i} = d{i} - vs{i}'*b{i};  %                             Step 2

v{i}(P,:) = U\(L\b{i});          % Solve A(:,P)*v(P,:) = b,    Step 3
delta{i} = d{i} - bs{i}'*v{i};   %                             Step 4

if m==2
  i = 2;
  
  M{i}  = B(1:end-m+i,1:end-m+i);
  b{i}  = M{i}(1:end-1,end);
  bs{i} = M{i}(end,1:end-1)';
  d{i}  = M{i}(end,end);

  v{i} = BEMsolve(i-1,b{i},v,vs,delta,deltas,b,bs,d,L,U,P,1);  % Step 1
  delta{i} = d{i} - v{i}'*bs{i};                               % Step 2 

  vs{i} = BEMsolve(i-1,bs{i},vs,v,deltas,delta,b,bs,d,L,U,P,2);% Step 3
  deltas{i} = d{i} - vs{i}'*b{i};                              % Step 4

end

x = BEMsolve(i,f,v,vs,delta,deltas,b,bs,d,L,U,P,1);
% JH: Causing Problems 10/30/06
% rel_resB = norm(B*x - r)/norm(r);  %  test
% if rel_resB > cis.LinSystResTol    %  test
%   %rel_resB;                         %  test ?
%   condB = condest(B);               %  test ?
%   if condB > cis.LinSystCondTol
%     fail = 1;
%   end
% end


%------------------------------------------------------

function x = BEMsolve(j,fi,v,vs,delta,deltas,b,bs,d,L,U,P,ot)

f{j} = fi;

for i=j:-1:1
  y{i} = [-vs{i}; 1]'*f{i}/deltas{i};  % Step 1
  fg   = f{i} - [b{i}; d{i}]*y{i};     % Step 2
  if i == 1
    f0     = fg(1:end-1,:);
  else
    f{i-1} = fg(1:end-1,:);
  end
  g{i}   = fg(end,:);
end

if ot == 1
  x0(P,:) =U\(L\f0);
else
  x0 = L'\(U'\f0(P,:));
end

for i=1:j
  if i == 1
    xt = x0;
  else
    xt = x{i-1};
  end
  ypp= (g{i} - bs{i}'*xt)/delta{i};     % Step 4
  x{i}  = [xt; y{i}] +  [-v{i}; 1]*ypp;
end

x = x{j};  
 
