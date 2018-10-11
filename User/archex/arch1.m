function bending

E = 1e3;
nu = 0.3;
boundmesh(-12, -12, 12, 12);

% Beam mesh
[nodex, ix, Iactive, tipy_id] = arch_mesh;
nodeu = 0*nodex;

% Load stepping
prev_tip = 0;
for tip = 0.1:0.1:24
  fprintf('-- Center = %g --\n', tip);

  % Predictor
  dtip = tip-prev_tip;
  [f1,K] = eval_Kf(nodeu, nodex, ix, E, nu);
  du = K(Iactive,Iactive)\(dtip * K(Iactive,tipy_id));
  nodeu(Iactive) = nodeu(Iactive) + du;
  nodeu(tipy_id) = -tip;
  prev_tip = tip;

  % Corrector
  for step = 1:5
    [f,K] = eval_Kf(nodeu, nodex, ix, E, nu);
    du = K(Iactive,Iactive)\f(Iactive);
    fprintf(' %d: %g %g\n', step, norm(f(Iactive)), norm(du));
    nodeu(Iactive) = nodeu(Iactive) + du;
  end

  K = (K+K')/2;
  [R,p] = chol(K(Iactive,Iactive));
  fprintf(' Stable: %d\n', (p==0));
  writemesh(nodex+nodeu, ix);
  pause(0.1);
end


function [nodex, ix, Iactive, tipy_id] = arch_mesh

[ix, nodex] = block2d9(pi,9.75, 0,10.75, 15,1);
idf     = reshape(1:prod(size(nodex)), size(nodex));
tipy_id = idf(2,48);
idf(:, 2) = 0;
idf(:,92) = 0;
idf(tipy_id) = 0;
Iactive = find(idf);
nodex = [cos(nodex(1,:)).*nodex(2,:);
	 sin(nodex(1,:)).*nodex(2,:)];

