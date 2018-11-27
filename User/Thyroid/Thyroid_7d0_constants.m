% Thyroid_5d_0_constants.m


%k <- c(v0=1,v1=1,v2=1,v3=1,a31=2.6e-2, a32=1.3e5, as=0.4, as2=2.6e5, at=0.1, b31=8e-6, b32=8.3e-4, bs=2.3e-4, bs2=140, bt=1.1e-6, DH=4.7e-8, DR=1e-10, Ds=50, DT=2.75, GD1=2.2e-8, GD2=4.3e-15, GT3=3.94e-13, GH=472, GR=1, GT=3.4e-12, km1=5e-7, km2=1e-9, Ls=1.68e6, Ss=100, TRH=6.9e-9, IBS=8e-6, TBG=3e-7, TBPA=4.5e-6, k=1, k30=2e9, k31=2e9, k41=2e10, k42=2e8);

%v0=1;
%v1=1;
%v01=1;
%v2=1;
%v3=1;
%v4=1; set v4 in front of TRH and vary v4 then from 0 to 10 
a31=2.6e-2; 
a32=1.3e5; 
as=0.4; 
as2=2.6e5; 
at=0.1; 
b31=8e-6; 
b32=8.3e-4; 
bs=2.3e-4; 
bs2=140; 
bt=1.1e-6; 
DH=4.7e-8; 
DR=1e-10; 
Ds=50; 
DT=2.75; 
GD1=2.2e-8; 
GD2=4.3e-15; 
GT3=3.94e-13; 
GH=472; 
GR=1; 
GT=3.4e-12; 
km1=5e-7; 
km2=1e-9; 
Ls=1.68e6; 
%Ss=100; 
%TRH0=6.9e-9
TRH0=6.9e-9; 
IBS=8e-6; 
TBG=3e-7; 
TBPA=4.5e-6; 
k=1; 
k30=2e9; 
k31=2e9; 
k41=2e10; 
k42=2e8;

% Original initial data Thyroid_5d_0_IVP initial values
%g = [1.7791e-11 5.2709e-12 1.1772e-08 1.8121e+00 1.9350e+00 0.36 0.36]';
y0_orig = [1.7791e-11 5.2709e-12 1.1772e-08 1.8121e+00 1.9350e+00 0 0]';
%
%optimum values (Johannes), normalizing factors: 
y1 = 17.7e-12;
y2 = 5.4e-12;
y3 = 1.2e-8;
y4 = 1.8;
y5 = 2.2;
y6 = 1;
y7 = 1;
ys = [y1 y2 y3 y4 y5 y6 y7]';
%
% scaled initial data
y0_scaled = y0_orig./ys;

