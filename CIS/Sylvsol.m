% Dense Sylvester equation solver

function [X,evl] = Sylvsol(A,B,C)

[X,evl] = lyap(A,-B,-C);     % Solve Sylvester equation by lyapa

%--< END OF Sylvsol >--


function [X,evl] = lyap(A, B, C)
%LYAP  Solve continuous-time Lyapunov equations.
%
%   X = LYAP(A,B,C) solves the general form of the Lyapunov matrix
%   equation (also called Sylvester equation):
%
%           A*X + X*B = -C
%
%   See also  DLYAP.

%        S.N. Bangert 1-10-86
%        Copyright 1986-2000 The MathWorks, Inc.
%        $Revision: 1.1 $  $Date: 2003/02/09 19:38:26 $
%        Last revised JNL 3-24-88, AFP 9-3-95

%ni = nargin;

[ma,na] = size(A);
[mb,nb] = size(B);
[mc,nc] = size(C);

% A and B must be square and C must have the rows of A and columns of B
if (ma ~= na) || (mb ~= nb) || (mc ~= ma) || (nc ~= mb)
   error('Dimensions do not agree.')
elseif ma==0 || mb==0,
   X = zeros(ma,mb); evl = zeros(ma+mb, 1); % DV
   return
end

% Perform schur decomposition on A (and convert to complex form)
[ua,ta] = schur(A);
[ua,ta] = rsf2csf(ua,ta);

% Perform schur decomposition on B (and convert to complex form)
[ub,tb] = schur(B);
[ub,tb] = rsf2csf(ub,tb);
evl = [diag(-tb); diag(ta)]; % get the eigenvalues

% Check all combinations of ta(i,i)+tb(j,j) for zero
p1 = diag(ta).'; % Use .' instead of ' in case A and B are not real
p1 = p1(ones(mb,1),:);
p2 = diag(tb);
p2 = p2(:,ones(ma,1));
sum = abs(p1) + abs(p2);
if any(any(sum == 0)) || any(any(abs(p1 + p2) < 1000*eps*sum))
   error('Solution does not exist or is not unique.');
end

% Transform C
ucu = -ua'*C*ub;

% Solve for first column of transformed solution
y = zeros(ma,mb);
ema = eye(ma);
y(:,1) = (ta+ema*tb(1,1))\ucu(:,1);

% Solve for remaining columns of transformed solution
for k=2:mb,
   km1 = 1:(k-1);
   y(:,k) = (ta+ema*tb(k,k))\(ucu(:,k)-y(:,km1)*tb(km1,k));
end

% Find untransformed solution
X = ua*y*ub';

% Ignore complex part if real inputs (better be small)
if isreal(A) && isreal(B) && isreal(C),
   X = real(X);
end

%--< END OF lyap >--
