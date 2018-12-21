[x75,v75,~,mult75] = loadPoint(...
  fullfile('Data','fusion_Orb_LC_N_75_16-Dec-2018_16_58_59.dat'));
x0 = x75(:,416);
v0 = v75(:,416);

N = 75;                     
odefile = str2func(sprintf('fusion_precomputed_with_sage_N_%d', N));
a = -1;
b = -0.3;
q_inf = -0.72;
parameters = {a ; b; q_inf};

%nPhases = 3*(N-1);
period = x0(end-1);
q_inf = x0(end);
parameters = {a ; b; q_inf};

init_LC_LC_L(odefile, x0, v0, s, par, ap,ntst,ncol)