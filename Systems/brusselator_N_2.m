function out = brusselator_N_2
out{1} = @init;
out{2} = @fun_eval;
out{3} = @jacobian;
out{4} = @jacobianp;
out{5} = @hessians;
out{6} = @hessiansp;
out{7} = [];
out{8} = [];
out{9} = [];
out{10} = [];

% --------------------------------------------------------------------------
function dydt = fun_eval(t,kmrgd,par_L,par_A,par_B,par_Dx,par_Dy)
dydt=[par_A-kmrgd(1)*(par_B+1)+kmrgd(1)^2*kmrgd(3)+(9*par_Dx*(par_A-2*kmrgd(1)+kmrgd(2)))/par_L^2;;
par_A-kmrgd(2)*(par_B+1)+kmrgd(2)^2*kmrgd(4)+(9*par_Dx*(par_A+kmrgd(1)-2*kmrgd(2)))/par_L^2;;
par_B*kmrgd(1)-kmrgd(1)^2*kmrgd(3)+(9*par_Dy*(kmrgd(4)-2*kmrgd(3)+par_B/par_A))/par_L^2;;
par_B*kmrgd(2)-kmrgd(2)^2*kmrgd(4)+(9*par_Dy*(kmrgd(3)-2*kmrgd(4)+par_B/par_A))/par_L^2;;];

% --------------------------------------------------------------------------
function [tspan,y0,options] = init
handles = feval(brusselator_N_2);
y0=[0,0,0,0];
options = odeset('Jacobian',handles(3),'JacobianP',handles(4),'Hessians',handles(5),'HessiansP',handles(6));
tspan = [0 10];

% --------------------------------------------------------------------------
function jac = jacobian(t,kmrgd,par_L,par_A,par_B,par_Dx,par_Dy)
jac=[ 2*kmrgd(1)*kmrgd(3) - par_B - (18*par_Dx)/par_L^2 - 1 , (9*par_Dx)/par_L^2 , kmrgd(1)^2 , 0 ; (9*par_Dx)/par_L^2 , 2*kmrgd(2)*kmrgd(4) - par_B - (18*par_Dx)/par_L^2 - 1 , 0 , kmrgd(2)^2 ; par_B - 2*kmrgd(1)*kmrgd(3) , 0 , - (18*par_Dy)/par_L^2 - kmrgd(1)^2 , (9*par_Dy)/par_L^2 ; 0 , par_B - 2*kmrgd(2)*kmrgd(4) , (9*par_Dy)/par_L^2 , - (18*par_Dy)/par_L^2 - kmrgd(2)^2 ];
% --------------------------------------------------------------------------
function jacp = jacobianp(t,kmrgd,par_L,par_A,par_B,par_Dx,par_Dy)
jacp=[ -(18*par_Dx*(kmrgd(2) - 2*kmrgd(1) + par_A))/par_L^3 , (9*par_Dx)/par_L^2 + 1 , -kmrgd(1) , (9*(kmrgd(2) - 2*kmrgd(1) + par_A))/par_L^2 , 0 ; -(18*par_Dx*(kmrgd(1) - 2*kmrgd(2) + par_A))/par_L^3 , (9*par_Dx)/par_L^2 + 1 , -kmrgd(2) , (9*(kmrgd(1) - 2*kmrgd(2) + par_A))/par_L^2 , 0 ; -(18*par_Dy*(kmrgd(4) - 2*kmrgd(3) + par_B/par_A))/par_L^3 , -(9*par_B*par_Dy)/(par_A^2*par_L^2) , kmrgd(1) + (9*par_Dy)/(par_A*par_L^2) , 0 , (9*(kmrgd(4) - 2*kmrgd(3) + par_B/par_A))/par_L^2 ; -(18*par_Dy*(kmrgd(3) - 2*kmrgd(4) + par_B/par_A))/par_L^3 , -(9*par_B*par_Dy)/(par_A^2*par_L^2) , kmrgd(2) + (9*par_Dy)/(par_A*par_L^2) , 0 , (9*(kmrgd(3) - 2*kmrgd(4) + par_B/par_A))/par_L^2 ];
% --------------------------------------------------------------------------
function hess = hessians(t,kmrgd,par_L,par_A,par_B,par_Dx,par_Dy)
hess1=[ 2*kmrgd(3) , 0 , 2*kmrgd(1) , 0 ; 0 , 0 , 0 , 0 ; -2*kmrgd(3) , 0 , -2*kmrgd(1) , 0 ; 0 , 0 , 0 , 0 ];
hess2=[ 0 , 0 , 0 , 0 ; 0 , 2*kmrgd(4) , 0 , 2*kmrgd(2) ; 0 , 0 , 0 , 0 ; 0 , -2*kmrgd(4) , 0 , -2*kmrgd(2) ];
hess3=[ 2*kmrgd(1) , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ; -2*kmrgd(1) , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ];
hess4=[ 0 , 0 , 0 , 0 ; 0 , 2*kmrgd(2) , 0 , 0 ; 0 , 0 , 0 , 0 ; 0 , -2*kmrgd(2) , 0 , 0 ];
hess(:,:,1) =hess1;
hess(:,:,2) =hess2;
hess(:,:,3) =hess3;
hess(:,:,4) =hess4;
% --------------------------------------------------------------------------
function hessp = hessiansp(t,kmrgd,par_L,par_A,par_B,par_Dx,par_Dy)
hessp1=[ (36*par_Dx)/par_L^3 , -(18*par_Dx)/par_L^3 , 0 , 0 ; -(18*par_Dx)/par_L^3 , (36*par_Dx)/par_L^3 , 0 , 0 ; 0 , 0 , (36*par_Dy)/par_L^3 , -(18*par_Dy)/par_L^3 ; 0 , 0 , -(18*par_Dy)/par_L^3 , (36*par_Dy)/par_L^3 ];
hessp2=[ 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ];
hessp3=[ -1 , 0 , 0 , 0 ; 0 , -1 , 0 , 0 ; 1 , 0 , 0 , 0 ; 0 , 1 , 0 , 0 ];
hessp4=[ -18/par_L^2 , 9/par_L^2 , 0 , 0 ; 9/par_L^2 , -18/par_L^2 , 0 , 0 ; 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ];
hessp5=[ 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ; 0 , 0 , -18/par_L^2 , 9/par_L^2 ; 0 , 0 , 9/par_L^2 , -18/par_L^2 ];
hessp(:,:,1) =hessp1;
hessp(:,:,2) =hessp2;
hessp(:,:,3) =hessp3;
hessp(:,:,4) =hessp4;
hessp(:,:,5) =hessp5;
%---------------------------------------------------------------------------
function tens3  = der3(t,kmrgd,par_L,par_A,par_B,par_Dx,par_Dy)
%---------------------------------------------------------------------------
function tens4  = der4(t,kmrgd,par_L,par_A,par_B,par_Dx,par_Dy)
%---------------------------------------------------------------------------
function tens5  = der5(t,kmrgd,par_L,par_A,par_B,par_Dx,par_Dy)
