
if ~ exist('x1','var')
  x1=loadPoint('na_tubular_reactor_orb_lc_11-Apr-2019_14_39_06.dat');
  x2=loadPoint('na_tubular_reactor_orb_lc_11-Apr-2019_14_50_08.dat');
end
fig = figure;
hold on
plot(x1(end,:),x1(end-1,:),'b');
plot(x2(end,:),x2(end-1,:),'b');

xlabel('D')
ylabel('period')

load('na_tubular_reactor_orb_lc_11-Apr-2019_14_39_06.mat')

for i=2:length(s)
  plot_singularity(s(i))
end

load('na_tubular_reactor_orb_lc_11-Apr-2019_14_50_08.mat')
 
for i=2:length(s)
  plot_singularity(s(i))
end

saveas(fig,'/home/carel/Documents/Master Thesis/reactor N:100 with bifs.svg');

xlim([0.1595 0.161])

saveas(fig,'/home/carel/Documents/Master Thesis/reactor N:100 detail.svg');

