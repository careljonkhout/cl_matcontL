function [borders, index1, index2] = find_new_borders_H(CISdata, x0, p, kapa)

global cds contopts

NSub = CISdata.NSub;
T    = CISdata.T;
Q1   = CISdata.Q;
A1   = CISdata.A;
Asub = T(1:NSub,1:NSub);

if nargin < 4
    [~, d ] = eig(Asub);
    d = diag(d);
    [~,ind] = min(abs(real(d)));
    omega = imag(d(ind(1)));
    kapa  = omega^2;
else
    omega = sqrt(kapa);
end

%% determine borders
% borders.v
[V,d] = eig(Asub);
d = diag(d);
idx1 = find(abs(d-1i*omega) == min(abs(d-1i*omega)));
if imag(d(idx1)) == 0
    error('Neutral saddle');
end
[Q,~,~] = qr([real(V(:,idx1)) imag(V(:,idx1))]);
borders.v = Q(:,1:2);

% borders.w
[V,d] = eig(Asub');
d = diag(d);
idx1 = find(abs(d+1i*omega) == min(abs(d+1i*omega)));
[Q,~,~] = qr([real(V(:,idx1)) imag(V(:,idx1))]);
borders.w = Q(:,1:2);

%% determine indices of lowest condition number
RED  = Asub*Asub+kapa*eye(NSub);
Bord  = [RED borders.w; borders.v' zeros(2)];
bunit = [zeros(NSub,2);eye(2)];
vext  = Bord\bunit;
wext  = Bord'\bunit;

Q_partial = Q1(:,1:NSub); %test FQ 08/03/2006

% DV: Compute true left eigenvector of A^2 + kapa*I 
RED2 = A1*A1+kapa*speye(length(A1));
Bord2 = [RED2' Q_partial*borders.v; borders.w'*Q_partial' zeros(2)];
bunit2 = [zeros(length(A1), 2); eye(2)];
wext2 = bordCIS1(Bord2, bunit2, 2);    
Qw1   = wext2(1:end-2,1);
Qw2   = wext2(1:end-2,2);

v1 = vext(1:NSub,1);
v2 = vext(1:NSub,2);
w1 = wext(1:NSub,1);
w2 = wext(1:NSub,2);

Qv1  = Q_partial*v1;
Qv2  = Q_partial*v2;
QCv1 = Q_partial*Asub*v1;
QCv2 = Q_partial*Asub*v2;
% Qw1  = Q1_partial*w1;
% Qw2  = Q1_partial*w2;

gx(1,:) = -Qw1'*hjacDirDerSquare(x0,p,Q_partial,Asub*v1) - Qw1'*A1*hjacDirDerSquare(x0,p,Q_partial,v1);
gx(2,:) = -Qw1'*hjacDirDerSquare(x0,p,Q_partial,Asub*v2) - Qw1'*A1*hjacDirDerSquare(x0,p,Q_partial,v2);
gx(3,:) = -Qw2'*hjacDirDerSquare(x0,p,Q_partial,Asub*v1) - Qw2'*A1*hjacDirDerSquare(x0,p,Q_partial,v1);
gx(4,:) = -Qw2'*hjacDirDerSquare(x0,p,Q_partial,Asub*v2) - Qw2'*A1*hjacDirDerSquare(x0,p,Q_partial,v2);

gk(1,1) = -w1'*v1;
gk(2,1) = -w1'*v2;
gk(3,1) = -w2'*v1;
gk(4,1) = -w2'*v2;

nap = cds.nap;
ap  = cds.ActiveParams;
Incr = contopts.Increment;

gp = zeros(4, nap);
p_temp = cell2mat(p);
for i=1:nap
    p1 = p_temp;
    p1(ap(i)) = p1(ap(i)) - Incr;
    p1 = num2cell(p1);
    f_x1= hjac(x0,p1);
    p2 = p_temp;
    p2(ap(i)) = p2(ap(i)) + Incr;
    p2 = num2cell(p2);
    f_x2= hjac(x0,p2);

    gp(1,i) = -Qw1'*(f_x2*QCv1 - f_x1*QCv1)/(2*Incr) - Qw1'*A1*(f_x2*Qv1 - f_x1*Qv1)/(2*Incr);
    gp(2,i) = -Qw1'*(f_x2*QCv2 - f_x1*QCv2)/(2*Incr) - Qw1'*A1*(f_x2*Qv2 - f_x1*Qv2)/(2*Incr);
    gp(3,i) = -Qw2'*(f_x2*QCv1 - f_x1*QCv1)/(2*Incr) - Qw2'*A1*(f_x2*Qv1 - f_x1*Qv1)/(2*Incr);
    gp(4,i) = -Qw2'*(f_x2*QCv2 - f_x1*QCv2)/(2*Incr) - Qw2'*A1*(f_x2*Qv2 - f_x1*Qv2)/(2*Incr);
end

two_index = zeros(6,2);
two_index(1,:) = [1 2];
two_index(2,:) = [1 3];
two_index(3,:) = [1 4];
two_index(4,:) = [2 3];
two_index(5,:) = [2 4];
two_index(6,:) = [3 4];

transferForm = [1,1;1,2;2,1;2,2];

condition_number_partial = zeros(6,1);
jac  = hjac(x0,p);
jacp = hjacp(x0,p);
A    = [jac  jacp zeros(size(jac,1),1)];
matrix_length = size(A,2);
for index=1:6
    number1 = two_index(index,1);
    number2 = two_index(index,2);

    testMatrix = [A;gx(number1,:) gp(number1,:) gk(number1,:);gx(number2,:) gp(number2,:) gk(number2,:);eye(1,matrix_length)];
    condition_number_partial(index) = partial_cond(testMatrix,3);
end

[Y, I] = min(condition_number_partial);

Index1 = two_index(I(1),1);
Index2 = two_index(I(1),2);

index1 = transferForm(Index1,:);
index2 = transferForm(Index2,:);