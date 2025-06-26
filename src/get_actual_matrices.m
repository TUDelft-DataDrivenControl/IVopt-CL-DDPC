function [Up2Yf_inno, Yp2Yf_inno, Uf2Yf_inno] = get_actual_matrices(A,B,C,D,K,p,f)
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
% tT_u= make_blk_tril_toeplitz(At,Bt,C,D,pf_max);
% tTf_u = tT_u(1:f*ny,1:f*nu); 
% tTp_u = tT_u(1:p*ny,1:p*nu); clear tT_u;
% tTf_y = make_blk_tril_toeplitz(At,K,C,zeros(ny,ny),f);
% T_e = make_blk_tril_toeplitz(A,K,C,eye(ny,ny),pf_max);
% Tp_e = T_e(1:p*ny,1:p*ny);
% Tf_e = T_e(1:f*ny,1:f*ny); clear T_e;
% % note: inv(I - tTf_y) * [tGf tTf_u eye(ny*f)] == [Gf Tf_u Tf_e]
% % (eye(ny*f)-tTf_y)\[tGf tTf_u eye(ny*f)]-[Gf Tf_u Tf_e] % <- check ==0

% useful constants
Atp_pinvGp = At^p*pinv(Gp);
Up2X0 = tKp_u-Atp_pinvGp*Tp_u; % Up -> X0
Yp2X0 = tKp_y+Atp_pinvGp;      % Yp -> X0
% Ep2X0 = -Atp_pinvGp*Tp_e;      % Ep -> X0

%{
% predictor form matrices:
Up2Yf_pred = tGf*Up2X0;
Uf2Yf_pred = tTf_u;
Yp2Yf_pred = tGf*Yp2X0;
Yf2Yf_pred = tTf_y;
Ep2Yf_pred = tGf*Ep2X0;
Ef2Yf_pred = eye(ny*f,ny*f);
%}

% innovation form matrices:
Up2Yf_inno = Gf*Up2X0;
Uf2Yf_inno = Tf_u;
Yp2Yf_inno = Gf*Yp2X0;
% Ep2Yf_inno = Gf*Ep2X0;
% Ef2Yf_inno = Tf_e;
end