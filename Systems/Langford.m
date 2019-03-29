function out = Langford
out{1} = @init;
out{2} = @fun_eval;
out{3} = @jacobian;
out{4} = @jacobianp;
out{5} = @hessians;
out{6} = @hessiansp;
out{7} = @der3;
out{8} = [];
out{9} = [];
out{10} = [];

% --------------------------------------------------------------------------
function dydt = fun_eval(t,kmrgd,par_b,par_c,par_d,par_LAMBDA)
dydt=[(par_LAMBDA-par_b)*kmrgd(1)-par_c*kmrgd(2)+kmrgd(1)*kmrgd(3)+par_d*kmrgd(1)*(1-kmrgd(3)^2);
par_c*kmrgd(1)+(par_LAMBDA-par_b)*kmrgd(2)+kmrgd(2)*kmrgd(3)+par_d*kmrgd(2)*(1-kmrgd(3)^2);
par_LAMBDA*kmrgd(3)-(kmrgd(1)^2+kmrgd(2)^2+kmrgd(3)^2);];

% --------------------------------------------------------------------------
function [tspan,y0,options] = init
handles = feval(Langford);
y0=[0,0,0];
options = odeset('Jacobian',handles(3),'JacobianP',handles(4),'Hessians',handles(5),'HessiansP',handles(6));
tspan = [0 10];

% --------------------------------------------------------------------------
function jac = jacobian(t,kmrgd,par_b,par_c,par_d,par_LAMBDA)
jac=[ kmrgd(3) - par_b + par_LAMBDA - par_d*(kmrgd(3)^2 - 1) , -par_c , kmrgd(1) - 2*kmrgd(1)*kmrgd(3)*par_d ; par_c , kmrgd(3) - par_b + par_LAMBDA - par_d*(kmrgd(3)^2 - 1) , kmrgd(2) - 2*kmrgd(2)*kmrgd(3)*par_d ; -2*kmrgd(1) , -2*kmrgd(2) , par_LAMBDA - 2*kmrgd(3) ];
% --------------------------------------------------------------------------
function jacp = jacobianp(t,kmrgd,par_b,par_c,par_d,par_LAMBDA)
jacp=[ -kmrgd(1) , -kmrgd(2) , -kmrgd(1)*(kmrgd(3)^2 - 1) , kmrgd(1) ; -kmrgd(2) , kmrgd(1) , -kmrgd(2)*(kmrgd(3)^2 - 1) , kmrgd(2) ; 0 , 0 , 0 , kmrgd(3) ];
% --------------------------------------------------------------------------
function hess = hessians(t,kmrgd,par_b,par_c,par_d,par_LAMBDA)
hess1=[ 0 , 0 , 1 - 2*kmrgd(3)*par_d ; 0 , 0 , 0 ; -2 , 0 , 0 ];
hess2=[ 0 , 0 , 0 ; 0 , 0 , 1 - 2*kmrgd(3)*par_d ; 0 , -2 , 0 ];
hess3=[ 1 - 2*kmrgd(3)*par_d , 0 , -2*kmrgd(1)*par_d ; 0 , 1 - 2*kmrgd(3)*par_d , -2*kmrgd(2)*par_d ; 0 , 0 , -2 ];
hess(:,:,1) =hess1;
hess(:,:,2) =hess2;
hess(:,:,3) =hess3;
% --------------------------------------------------------------------------
function hessp = hessiansp(t,kmrgd,par_b,par_c,par_d,par_LAMBDA)
hessp1=[ -1 , 0 , 0 ; 0 , -1 , 0 ; 0 , 0 , 0 ];
hessp2=[ 0 , -1 , 0 ; 1 , 0 , 0 ; 0 , 0 , 0 ];
hessp3=[ 1 - kmrgd(3)^2 , 0 , -2*kmrgd(1)*kmrgd(3) ; 0 , 1 - kmrgd(3)^2 , -2*kmrgd(2)*kmrgd(3) ; 0 , 0 , 0 ];
hessp4=[ 1 , 0 , 0 ; 0 , 1 , 0 ; 0 , 0 , 1 ];
hessp(:,:,1) =hessp1;
hessp(:,:,2) =hessp2;
hessp(:,:,3) =hessp3;
hessp(:,:,4) =hessp4;
%---------------------------------------------------------------------------
function tens3  = der3(t,kmrgd,par_b,par_c,par_d,par_LAMBDA)
tens31=[ 0 , 0 , 0 ; 0 , 0 , 0 ; 0 , 0 , 0 ];
tens32=[ 0 , 0 , 0 ; 0 , 0 , 0 ; 0 , 0 , 0 ];
tens33=[ 0 , 0 , -2*par_d ; 0 , 0 , 0 ; 0 , 0 , 0 ];
tens34=[ 0 , 0 , 0 ; 0 , 0 , 0 ; 0 , 0 , 0 ];
tens35=[ 0 , 0 , 0 ; 0 , 0 , 0 ; 0 , 0 , 0 ];
tens36=[ 0 , 0 , 0 ; 0 , 0 , -2*par_d ; 0 , 0 , 0 ];
tens37=[ 0 , 0 , -2*par_d ; 0 , 0 , 0 ; 0 , 0 , 0 ];
tens38=[ 0 , 0 , 0 ; 0 , 0 , -2*par_d ; 0 , 0 , 0 ];
tens39=[ -2*par_d , 0 , 0 ; 0 , -2*par_d , 0 ; 0 , 0 , 0 ];
tens3(:,:,1,1) =tens31;
tens3(:,:,1,2) =tens32;
tens3(:,:,1,3) =tens33;
tens3(:,:,2,1) =tens34;
tens3(:,:,2,2) =tens35;
tens3(:,:,2,3) =tens36;
tens3(:,:,3,1) =tens37;
tens3(:,:,3,2) =tens38;
tens3(:,:,3,3) =tens39;
%---------------------------------------------------------------------------
function tens4  = der4(t,kmrgd,par_b,par_c,par_d,par_LAMBDA)
%---------------------------------------------------------------------------
function tens5  = der5(t,kmrgd,par_b,par_c,par_d,par_LAMBDA)
