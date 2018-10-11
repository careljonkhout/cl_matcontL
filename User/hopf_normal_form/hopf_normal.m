function out = hopf_normal
out{1} = @init;
out{2} = @fun_eval;
out{3} = @jacobian;
out{4} = @jacobianp;
out{5} = @hessians;
out{6} = @hessiansp;
out{7} = @der3;
out{8} = @der4;
out{9} = @der5;
out{10} = [];

% --------------------------------------------------------------------------
function dydt = fun_eval(t,kmrgd,par_beta)
dydt = zeroes(20,1);
dydt(1)=par_beta*kmrgd(1)-kmrgd(2)+kmrgd(1)*(kmrgd(1)^2+kmrgd(2)^2);
dydt(2)=kmrgd(1)+par_beta*kmrgd(2)+kmrgd(2)*(kmrgd(1)^2+kmrgd(2)^2);];

% --------------------------------------------------------------------------
function [tspan,y0,options] = init
handles = feval(hopfnormal);
y0=[0,0];
options = odeset('Jacobian',handles(3),'JacobianP',handles(4),'Hessians',handles(5),'HessiansP',handles(6));
tspan = [0 10];

% --------------------------------------------------------------------------
function jac = jacobian(t,kmrgd,par_beta)
jac=[ par_beta + 3*kmrgd(1)^2 + kmrgd(2)^2 , 2*kmrgd(1)*kmrgd(2) - 1 ; 2*kmrgd(1)*kmrgd(2) + 1 , par_beta + kmrgd(1)^2 + 3*kmrgd(2)^2 ];
% --------------------------------------------------------------------------
function jacp = jacobianp(t,kmrgd,par_beta)
jacp=[ kmrgd(1) ; kmrgd(2) ];
% --------------------------------------------------------------------------
function hess = hessians(t,kmrgd,par_beta)
hess1=[ 6*kmrgd(1) , 2*kmrgd(2) ; 2*kmrgd(2) , 2*kmrgd(1) ];
hess2=[ 2*kmrgd(2) , 2*kmrgd(1) ; 2*kmrgd(1) , 6*kmrgd(2) ];
hess(:,:,1) =hess1;
hess(:,:,2) =hess2;
% --------------------------------------------------------------------------
function hessp = hessiansp(t,kmrgd,par_beta)
hessp1=[ 1 , 0 ; 0 , 1 ];
hessp(:,:,1) =hessp1;
%---------------------------------------------------------------------------
function tens3  = der3(t,kmrgd,par_beta)
tens31=[ 6 , 0 ; 0 , 2 ];
tens32=[ 0 , 2 ; 2 , 0 ];
tens33=[ 0 , 2 ; 2 , 0 ];
tens34=[ 2 , 0 ; 0 , 6 ];
tens3(:,:,1,1) =tens31;
tens3(:,:,1,2) =tens32;
tens3(:,:,2,1) =tens33;
tens3(:,:,2,2) =tens34;
%---------------------------------------------------------------------------
function tens4  = der4(t,kmrgd,par_beta)
tens41=[ 0 , 0 ; 0 , 0 ];
tens42=[ 0 , 0 ; 0 , 0 ];
tens43=[ 0 , 0 ; 0 , 0 ];
tens44=[ 0 , 0 ; 0 , 0 ];
tens45=[ 0 , 0 ; 0 , 0 ];
tens46=[ 0 , 0 ; 0 , 0 ];
tens47=[ 0 , 0 ; 0 , 0 ];
tens48=[ 0 , 0 ; 0 , 0 ];
tens4(:,:,1,1,1) =tens41;
tens4(:,:,1,1,2) =tens42;
tens4(:,:,1,2,1) =tens43;
tens4(:,:,1,2,2) =tens44;
tens4(:,:,2,1,1) =tens45;
tens4(:,:,2,1,2) =tens46;
tens4(:,:,2,2,1) =tens47;
tens4(:,:,2,2,2) =tens48;
%---------------------------------------------------------------------------
function tens5  = der5(t,kmrgd,par_beta)
tens51=[ 0 , 0 ; 0 , 0 ];
tens52=[ 0 , 0 ; 0 , 0 ];
tens53=[ 0 , 0 ; 0 , 0 ];
tens54=[ 0 , 0 ; 0 , 0 ];
tens55=[ 0 , 0 ; 0 , 0 ];
tens56=[ 0 , 0 ; 0 , 0 ];
tens57=[ 0 , 0 ; 0 , 0 ];
tens58=[ 0 , 0 ; 0 , 0 ];
tens59=[ 0 , 0 ; 0 , 0 ];
tens510=[ 0 , 0 ; 0 , 0 ];
tens511=[ 0 , 0 ; 0 , 0 ];
tens512=[ 0 , 0 ; 0 , 0 ];
tens513=[ 0 , 0 ; 0 , 0 ];
tens514=[ 0 , 0 ; 0 , 0 ];
tens515=[ 0 , 0 ; 0 , 0 ];
tens516=[ 0 , 0 ; 0 , 0 ];
tens5(:,:,1,1,1,1) =tens51;
tens5(:,:,1,1,1,2) =tens52;
tens5(:,:,1,1,2,1) =tens53;
tens5(:,:,1,1,2,2) =tens54;
tens5(:,:,1,2,1,1) =tens55;
tens5(:,:,1,2,1,2) =tens56;
tens5(:,:,1,2,2,1) =tens57;
tens5(:,:,1,2,2,2) =tens58;
tens5(:,:,2,1,1,1) =tens59;
tens5(:,:,2,1,1,2) =tens510;
tens5(:,:,2,1,2,1) =tens511;
tens5(:,:,2,1,2,2) =tens512;
tens5(:,:,2,2,1,1) =tens513;
tens5(:,:,2,2,1,2) =tens514;
tens5(:,:,2,2,2,1) =tens515;
tens5(:,:,2,2,2,2) =tens516;
function usernorm

