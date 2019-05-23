
close all
load('arguments_for_multipliers.mat')
global lds

JJ = [J;rand(1,(lds.ntst*lds.ncol+1)*lds.nphase + 2)];
b = rand((lds.ntst*lds.ncol+1)*lds.nphase+2,1);

tic; x = linear_solver_collocation_2(JJ,b); toc
tic; x_std = JJ \ b; toc

hold on
abs_diff = abs(x-x_std);
x = abs(x);
x_std = abs(x_std);


for i=1:2*lds.nphase
  close all
  figure
  hold on
coords = i:lds.ncol*lds.nphase:length(x)-2;

plot(abs_diff(coords));
plot(x    (coords));
plot(x_std(coords));
legend
pause
end



bb = JJ*x;

