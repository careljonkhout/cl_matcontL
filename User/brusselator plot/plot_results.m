
filename_lc_lc_1 = 'bruss_ep_h_1';
filename = filename_lc_lc_1;
datafile = fullfile([filename '.dat']);
x = loadPoint(datafile);


figure
hold on
plot(x(end,:),x(end-1,:))
xlabel('L')
ylabel('period')

load('s_bruss_oc')

ns_1_L = s(2).data.parametervalues(2);
ns_1_T = s(2).data.T;

plot(ns_1_L,ns_1_T,'r*')
text(ns_1_L,ns_1_T,'NS','VerticalAlignment','top')

ns_2_L = s(3).data.parametervalues(2);
ns_2_T = s(3).data.T;

plot(ns_2_L,ns_2_T,'r*')
text(ns_2_L,ns_2_T,'NS','VerticalAlignment','top')

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

