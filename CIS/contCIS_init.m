function CISdata = contCIS_init(A1, i_PT, NSub, NUnstable)
print_diag(5,'In contCIS_init\n');
% Reinitialization of CIS

global cds

if i_PT > 0
  NUnstable = -1;   % Recompute NUnstable
  cds.tfUpdate  =  1;   % need update test functions in MATCONT
  print_diag(3,'contCIS_init: need adapt test functions\n')           % CIS test
end

%JH: Modified second argument from constant 1 
[Q1,T1,evl1,NSub,NUnstable] = CISinit(A1, NSub, NUnstable);

[evl1_r, evl1_l] = SortEvl(evl1, NSub);

%DV: create CIS data struct
CISdata.A = A1;
CISdata.Q = Q1;
CISdata.T = T1;
CISdata.evl_r = evl1_r;
CISdata.evl_l = evl1_l;
CISdata.NSub = NSub;
CISdata.NUnstable = NUnstable;
CISdata.overlap = 0;