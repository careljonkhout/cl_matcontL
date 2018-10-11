function [f,K] = eval_Kf(nodeu, nodex, ix, E, nu)

nshape = size(ix,1);
nelt = size(ix,2);
nshape2 = nshape*2;
Ki = zeros(nshape2*nshape2, nelt);
Kj = zeros(nshape2*nshape2, nelt);
Ka = zeros(nshape2*nshape2, nelt);

N = prod(size(nodex));
f = zeros(N,1);
idf = reshape(1:N,size(nodex));
e   = ones(nshape2,1);

for ielt = 1:nelt
  
  % -- Get element contribution --
  I = ix(:,ielt);
  enodex = nodex(:,I);
  enodeu = nodeu(:,I);
        mex_id_ = 'element_neohook(i double[xx], i double[xx], i double, i double, o double[xx], o double[x], i int)';
[Kelt, felt] = femex(mex_id_, enodex, enodeu, E, nu, nshape, 2, nshape, 2, nshape, nshape2, nshape2, nshape2);
  
  % -- Assemble into K and f
  J = reshape(idf(:,I), nshape2, 1);
  f(J) = f(J) + felt;
  Ki(:,ielt) = reshape(e*J', nshape2*nshape2, 1);
  Kj(:,ielt) = reshape(J*e', nshape2*nshape2, 1);
  Ka(:,ielt) = reshape(Kelt, nshape2*nshape2, 1);
  
end

K = sparse(Ki, Kj, Ka, N, N);
