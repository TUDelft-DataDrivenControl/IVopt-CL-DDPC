function [Up2Yf_inno, Yp2Yf_inno, Uf2Yf_inno] = get_actual_matrices(plant,p,f)
%% Computes true components of the matrix Lf from known system dynamics
% Lf = [Up2Yf_inno Yp2Yf_inno Uf2Yf_inno]

[A,B,C,D,K] = plant2ABCDK(plant);
[ny,nu] = size(D);

At = A-K*C;
Bt = B-K*D;

pf_max = max(p,f);

% Yf =  Gf * X0 +  Tf_u * Uf +              Tf_e * Ef
% Yf = tGf * X0 + tTf_u * Uf + tTf_y * Yf +        Ef

% make extended observability matrices
G = make_ext_obsv(A,C,pf_max);
Gp = G(1:p*ny,:);
Gf = G(1:f*ny,:); clear G;
% tGf = make_ext_obsv(At,C,f);

% make extended reversed controllability matrices
tKp_u = make_ext_ctrb(At,Bt,p,rev=true);
tKp_y = make_ext_ctrb(At,K,p,rev=true);

% make block-lower toeplitz matrices
T_u = make_blk_tril_toeplitz(A,B,C,D,pf_max);
Tf_u = T_u(1:f*ny,1:f*nu); 
Tp_u = T_u(1:p*ny,1:p*nu); clear T_u;

% useful constants
Atp_pinvGp = At^p*pinv(Gp);
Up2X0 = tKp_u-Atp_pinvGp*Tp_u; % Up -> X0
Yp2X0 = tKp_y+Atp_pinvGp;      % Yp -> X0

% innovation form matrices:
Up2Yf_inno = Gf*Up2X0;
Uf2Yf_inno = Tf_u;
Yp2Yf_inno = Gf*Yp2X0;

end