
if ~ exist('x1','var')
  x1=loadPoint('lpc_fusion_oc_11-Apr-2019_12_12_55.dat');
  x2=loadPoint('lpc_fusion_oc_11-Apr-2019_12_47_40.dat');
end
fig = figure;
hold on
plot(x1(end,:),x1(end-1,:),'b');
plot(x2(end,:),x2(end-1,:),'b');

xlabel('q_{inf}')
ylabel('period')

load('/home/carel/Documents/cl_matcontL/User/fusion plot/lpc_fusion_oc_11-Apr-2019_14_05_42.mat')

for i=2:6
  plot_singularity(s(i))
end

load('/home/carel/Documents/cl_matcontL/User/fusion plot/lpc_fusion_oc_11-Apr-2019_12_47_40.mat')
 
for i=2:6
  plot_singularity(s(i))
end

saveas(fig,'/home/carel/Documents/Master Thesis/fusion N:25 with bifs.svg');
