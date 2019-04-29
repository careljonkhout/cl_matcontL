function set_ntst_ncol(ntst,ncol,newmsh)
% todo newm
%
% This function sets the number of mesh and collocation points
% as well as a new mesh distribution. All things that depend on
% the number of mesh and/or collocation points are updated.
%
global cds lds
lds.ntst = ntst;
lds.ncol = ncol;

lds.tsts = 1:ntst;
lds.cols = 1:ncol;
lds.tps = lds.ntst*lds.ncol+1; % number of points on curve

lds.ncoords = lds.tps*lds.nphase;
lds.coords = 1:lds.ncoords;
lds.PeriodIdx = lds.ncoords+1;
lds.idxmat = reshape(fix((1:((lds.ncol+1)*lds.ntst))/(1+1/lds.ncol))+1,lds.ncol+1,lds.ntst);
cds.ndim = lds.ncoords+2;
cds.oldJacX = [];

lds.msh = newmsh;
lds.finemsh = zeros(1,ntst*ncol+1);
fine_mesh_index = 2;

for i=1:ntst
  mesh_interval_width = lds.msh(i+1) - lds.msh(i);
  fine_mesh_dt        = mesh_interval_width / ncol;
  for j=1:ncol
    lds.finemsh(fine_mesh_index) = lds.msh(i) + j * fine_mesh_dt;
    fine_mesh_index              = fine_mesh_index + 1;
  end
end

lds.dt = lds.msh(lds.tsts+1)-lds.msh(lds.tsts);

lds.upoldp = [];
lds.multipliersX = [];
lds.CalcMultipliers =0;
calc_weigths;
