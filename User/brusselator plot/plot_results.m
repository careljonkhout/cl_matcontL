
filename_lc_lc_1 = 'bruss_ep_h_1';
filename = filename_lc_lc_1;
datafile = fullfile([filename '.dat']);
x = loadPoint(datafile);


fig = figure
hold on
plot(x(end,:),x(end-1,:))
xlabel('L')
ylabel('period')

load('bruss_oc_orb_lc')

plot_singularity(s(2), 'VerticalAlignment', 'bottom');
plot_singularity(s(3), 'VerticalAlignment', 'bottom');
plot_singularity(s(4), 'VerticalAlignment', 'bottom');
plot_singularity(s(5), 'VerticalAlignment', 'bottom');


filename_lc_lc_1 = 'bruss_ep_h_2';
filename = filename_lc_lc_1;
datafile = fullfile([filename '.dat']);
x = loadPoint(datafile);

plot(x(end,:),x(end-1,:));

load('bruss_ep_h_2')

bpc_L = s(2).data.parametervalues(2);
bpc_T = s(2).data.T;

plot(bpc_L,bpc_T,'r*')
text(bpc_L,bpc_T,'BPC','VerticalAlignment','top')

datafile = fullfile('bruss_bpc_lc.dat');
x = loadPoint(datafile);
plot(x(end,:),x(end-1,:))

load('bruss_bpc_lc.mat')
struct2table(s)
plot_singularity(s(4),'VerticalAlignment','top');
plot_singularity(s(5),'VerticalAlignment','top');



bpc_point = x(:,1);
L_bpc = bpc_point(end);
T_bpc = bpc_point(end-1);

plot(L_bpc,T_bpc,'r*')
text(L_bpc,T_bpc,'BPC','VerticalAlignment','bottom')


saveas(fig,'/home/carel/Documents/Master Thesis/bruss N:31.svg');

function plot_singularity(s, varargin)
  L = s.data.parametervalues(2);
  T = s.data.T;
  plot(L,T,'r*')
  text(L,T,s.label,varargin{:})
end





