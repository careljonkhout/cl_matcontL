function borders = find_newborders_LP(CISdata)

NSub   = CISdata.NSub;
Asub   = CISdata.T(1:NSub,1:NSub);
d = eig(Asub);
% find a real eigenvalue with smallest absolute value and call it Rmin
I    = find(imag(d) == 0);
Rmin = min(abs((d(I))));
I1   = find(abs((d)) == Rmin);
RED  = Asub - d(I1(1))*eye(NSub);
borders.v = null(RED);
borders.w = null(RED');

Bord=[Asub borders.w;borders.v' 0];
bunit=[zeros(NSub,1); 1];
vext=Bord\bunit;
wext=Bord'\bunit;
borders.v = vext(1:NSub)/norm(vext(1:NSub));
borders.w = wext(1:NSub)/norm(wext(1:NSub));