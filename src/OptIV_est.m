%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           DDPC using an Optimal-IV
%           Authors: R. Dinkla, T. Oomen, J.W. van Wingerden
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solves a constrained least-squares regression problem
clear;
close all;
rng default;

addpath('..\bin');

% color maps: https://www.fabiocrameri.ch/colourmaps/

%% problem parameters
f = 10;
fid = f;
p = 10;
ny = 1;
nu = 1;
N = 1e4;
Nbar = p + fid + N -1;

% controller weights
dRk= 1*eye(nu);  dR = kron(speye(f),dRk);
Rk = 1*eye(nu);   R = kron(speye(f),Rk);
Qk = 1e3*eye(ny); Q = kron(speye(f),Qk);

% reference
y_ref = idinput(Nbar-fid+2*f-1,'prbs',[],[-1 1]).';
% y_ref = sign(sin((1:3:3*(Nbar-fid+2*f))*2*pi/(Nbar/20)));
% y_ref = y_ref(:,1:end-1);
y_ref = repmat(y_ref,ny,1);

% get solver: uf_ref, yf_ref, up, yp -> u_k
get_solver;

%% ============================= make system ==============================
% make nominal system
% nx = 2;
% sys = drss(nx,ny,nu);% sys = c2d(sys,1e-1);
% [A,B,C,D] = ssdata(sys); D = 0*D; sys = ss(A,B,C,D,-1);
% model_Favoreel1999; % D = 0...
model_DoubleTank; R_e = 1e0*Re;

% set noise properties
% R_e = 1e-2*tril(rand(ny)*2);
% if any(diag(R_e)==0) % ensure Re > 0
%     var_e_diag = diag(R_e);
%     var_e_diag(var_e_diag==0) = 1e-2;
%     indiag = 1:ny;
%     R_e(indiag + (indiag - 1)*ny) = var_e_diag;
% end
% R_e = R_e*R_e.';
% K = place(A.',C.',linspace(0.1,0.15,nx)).';

% create system with noise
sys2 = ss(A,[B K],C,[D eye(ny)],-1);

% ------------------- calculate optimal steady state ----------------------
% P0 = C*((eye(nx)-A)\B)+D;
P0 = dcgain(sys2); P0 = P0(:,1:nu);
u_ref = P0\y_ref;

%% ================= Run #0: initial simulation (OL) ======================
R_u = eye(nu)*1e2; % set OL-input variance
u = mvnrnd(zeros(nu,1),R_u,Nbar).';
e = mvnrnd(zeros(ny,1),R_e,Nbar).';
[y,~,x] = lsim(sys2,[u;e],[]); y = y.'; x = x.';

%% ==================== Create initial controller =========================
% create Hankel matrices
[Upf_r0,Up_r0,Uf_r0] = make_Hankel(u,p,fid);
[Ypf_r0,Yp_r0,Yf_r0] = make_Hankel(y,p,fid);

Lest_r0 = Yf_r0*pinv([Up_r0;Yp_r0;Uf_r0]);
% [~, ~, ~, ~, Gf_tKp_u_r0,Tf_u_r0,Gf_tKp_y_r0,~] = L1est2mats(L1est_r0,p,f,nu);

% make estimate causal: block-lower triangularize Tf_u (influence uf -> yf)
% Tf_u_r0 = mat2cell(Tf_u_r0,ny*ones(1,f),nu*ones(1,f));
% Tf_u_r0(logical(triu(ones(f,f),1))) = {zeros(ny,nu)};
% Tf_u_r0 = cell2mat(Tf_u_r0);

% calculate matrices for analytical solution
% rf2u = (R+Tf_u_r0.'*Q*Tf_u_r0)\Tf_u_r0.'*Q;
% rf2u = rf2u(1:nu,:);
% up2u = -rf2u*Gf_tKp_u_r0;
% yp2u = -rf2u*Gf_tKp_y_r0;
Cz0 = @(up,yp,yrf,urf) calc_u_v2(up,yp,Lest_r0,yrf,urf);

%% ================= Run #1: Closed-loop data collection ==================

% initialize closed-loop simulation
up1 = u(:,end-p+1:end);                      % up
yp1 = y(:,end-p+1:end);                      % yp
e1 = mvnrnd(zeros(ny,1),R_e,Nbar).';         % noise
x1_0 = A*x(:,end) + B*u(:,end) + K*e(:,end); % initial state

% run closed-loop simulation with real system
[u1,x1,y1] = run_real_cl(up1,yp1,x1_0,y_ref,u_ref,Cz0,  e1,ny,nu,nx,Nbar,f,A,B,C,D,K);
[u2,x2,y2] = run_real_cl(up1,yp1,x1_0,y_ref,u_ref,Cz0,0*e1,ny,nu,nx,Nbar,f,A,B,C,D,K); % <- for optimal IV

% scaling data & system
uSF = diag(mean(abs(u1),2));
u1 = uSF\u1; u2 = uSF\u2; up1 = uSF\up1; u_ref = uSF\u_ref;
B = B*uSF; D = D*uSF;
sys2 = ss(A,[B K],C,[D eye(ny)],-1);

% creating Hankel matrices
[~,Up_r1,Uf_r1] = make_Hankel(u1,p,fid); [~,Up_r2,Uf_r2] = make_Hankel(u2,p,fid);
[~,Yp_r1,Yf_r1] = make_Hankel(y1,p,fid); [~,Yp_r2,Yf_r2] = make_Hankel(y2,p,fid); 
[~,Ep_r1,Ef_r1] = make_Hankel(e1,p,fid);
[~,yuRp,yuRf] = make_Hankel([y_ref;u_ref],p,2*f-1);
[~,yRp,yRf] = make_Hankel(y_ref,p,2*f-1);
[~,uRp,uRf] = make_Hankel(u_ref,p,2*f-1);
Wp_r1 = [Up_r1;Uf_r1;Yp_r1];

%% make actual matries and scaled controller components
up2uk = up2usol(Lest_r0);
yp2uk = uSF\yp2usol(Lest_r0);
yrf2uk = uSF\yrf2usol(Lest_r0);
urf2uk = urf2usol(Lest_r0);

% --------------------- create actual matrices ----------------------------
get_actual_matrices;

% ---------------- make transfer-function of controller -------------------
% Cz0yk  = tf([0 fliplr(yp2uk)],[1 -fliplr(up2uk)],1,'Variable','z^-1'); % y_k -> u_k
% Cz0yrf = tf(fliplr(yrf2uk),[1 -fliplr(up2uk)],1,'Variable','z^-1');    % yr_{k+f-1} -> u_k
% Cz0urf = tf(fliplr(urf2uk),[1 -fliplr(up2uk)],1,'Variable','z^-1');    % ur_{k+f-1} -> u_k
% 
% qinv = tf([0 1],1,-1,'Variable','q^-1');
% Mpnu = tf(ss(kron(ones(p,1),eye(nu))));
% Mpny = tf(ss(kron(ones(p,1),eye(ny))));
% Mfny = tf(ss(kron(ones(f,1),eye(ny))));
% Mfnu = tf(ss(kron(ones(f,1),eye(nu))));
% for k=0:p-1
%     Mpnu(1:end-nu*k,:) = qinv*Mpnu(1:end-nu*k,:);
%     Mpny(1:end-ny*k,:) = qinv*Mpny(1:end-ny*k,:);
% end
% for k=1:f-1
%     Mfny(1:end-ny*k,:) = qinv*Mfny(1:end-ny*k,:);
%     Mfnu(1:end-nu*k,:) = qinv*Mfnu(1:end-nu*k,:);
% end
% Rq = eye(nu)-up2uk*Mpnu;
% Sq = -yp2uk*Mpny;
% Tq = [yrf2uk*Mfny urf2uk*Mfnu]; % yr_{k+f-1} => u_k and ur_{k+f-1} => u_k

% --------------------- making state-space controller ---------------------
% x_k = [ u_[k-p,k-1]; y_[k-p,k-1]; yr_[k,k+f-2]; ur_[k,k+f-2] ]
% u_k = [ y_k; yr_{k+f-1}; ur_{k+f-1} ]
Cc = [up2uk yp2uk yrf2uk(:,1:end-ny) urf2uk(:,1:end-nu)];
Dc = [zeros(nu,ny) yrf2uk(:,end-ny+1:end) urf2uk(:,end-nu+1:end)];
Ic_up = eye((p-1)*nu); Ic_yp = eye((p-1)*ny);
Ac = blkdiag([zeros((p-1)*nu,nu) Ic_up],zeros(nu,ny),Ic_yp);
if f > 1
    Ic_yr = eye((f-2)*ny); Ic_ur = eye((f-2)*nu);
    Ac = blkdiag(Ac,zeros(ny),Ic_yr,zeros(ny,nu),Ic_ur);
    Ac = [Ac; zeros(nu,size(Ac,2))];
else
    Ac = [Ac; zeros(ny,p*(nu+ny))];
end
Ac((p-1)*nu+1:p*nu,:) = Cc;
Bc = [zeros((p-1)*nu,ny*2+nu);Dc;zeros((p-1)*ny,2*ny+nu)];
if f > 1
    Bc = [Bc;blkdiag(eye(ny),[zeros((f-2)*ny,ny);eye(ny)],[zeros((f-2)*nu,nu);eye(nu)])];
else
    Bc = [Bc; eye(ny) zeros(ny,nu+ny)];
end

Cz_ss = ss(Ac,Bc,Cc,Dc,-1);


%% plot data - w/ & w/o noise
close all;
fig1 = figure(1);
ax1_1=subplot(2,1,1);
plot(y2,'k','LineWidth',1.5); hold on; grid on;
plot(y_ref);
plot(y1,'r');

legend({'$y^{\mathrm{opt}}_{iv}$','ref','$y$'},'Interpreter','latex');
ax1_2=subplot(2,1,2);
plot(u2,'k'); hold on; grid on;
plot(u_ref);
plot(u1,'r');
legend({'$u^{\mathrm{opt}}_{iv}$','ref','$u$'},'Interpreter','latex');
linkaxes([ax1_1 ax1_2],'x');

% -------------------- Get system estimate using IV -----------------------
load('vik.mat');
num_iter = 2;
fro_EfUf= nan(1,num_iter+1);
fro_EfZ = nan(1,num_iter+1);
fro_EfUp = norm(Ef_r1*Up_r1.','fro')^2;
fro_EfYp = norm(Ef_r1*Yp_r1.','fro')^2;

% optimal IV matrix
Ziv_opt = [Up_r2;Uf_r2;Yp_r2];

% creating initial IV matrix
Ziv = [yuRp;yuRf;Yp_r1];

fro_EfZ(1) = norm(Ef_r1*Ziv.'/N,'fro')^2;

% plotting correlations and estimate of system parameters
fig10 = figure(10); 
Fig10FS = 12;
tl10 = tiledlayout(4+num_iter*1,4,'TileSpacing','compact','Padding','compact');

% Row 1  - plot real system parameters
nexttile;
nexttile;
nexttile;
imagesc_vik([Up2Yf_inno Uf2Yf_inno],vik,gca);
nexttile;
imagesc_vik(Yp2Yf_inno,vik,gca);

% Row 2
nexttile;
imagesc_parts(Ef_r1, vik, gca, yuRp, yuRf, Yp_r1);
title('$\frac{1}{N} E_f \mathcal{Z}^\top$', 'Interpreter', 'latex', 'FontSize', Fig10FS);
ylabel('$\mathcal{Z}=[R_p; R_f; Y_p]$', 'Interpreter', 'latex', 'FontSize', Fig10FS);
nexttile;
imagesc_parts(Wp_r1, vik, gca, yuRp, yuRf, Yp_r1);
title('$\frac{1}{N} W_p \mathcal{Z}^\top$, $W_p=[U_p; U_f; Y_p]$', 'Interpreter', 'latex', 'FontSize', Fig10FS);
nexttile; 
LestZr = Yf_r1 * Ziv.' / N * pinv(Wp_r1 * Ziv.' / N);
imagesc_vik(LestZr(:,1:end-ny*p), vik, gca);
title('$\frac{1}{N} Y_f\mathcal{Z}^\top (\frac{1}{N} W_p\mathcal{Z}^\top)^\dagger$, $U_{pf} \rightarrow Y_f$', 'Interpreter', 'latex', 'FontSize', Fig10FS);
nexttile;
imagesc_vik(LestZr(:,end-ny*p+1:end), vik, gca);
title('$\frac{1}{N} Y_f\mathcal{Z}^\top (\frac{1}{N} W_p\mathcal{Z}^\top)^\dagger$, $Y_p \rightarrow Y_f$', 'Interpreter', 'latex', 'FontSize', Fig10FS);

% Row 3
nexttile;
imagesc_parts(Ef_r1, vik, gca, Up_r1, Uf_r1, Yp_r1);
ylabel('$\mathcal{Z}=[U_p; U_f; Y_p]$', 'Interpreter', 'latex', 'FontSize', Fig10FS);
nexttile;
imagesc_parts(Wp_r1, vik, gca, Up_r1, Uf_r1, Yp_r1);
nexttile; 
LestWp = Yf_r1 * Wp_r1.' / N * pinv(Wp_r1 * Wp_r1.' / N);
imagesc_vik(LestWp(:,1:end-ny*p), vik, gca);
nexttile;
imagesc_vik(LestWp(:,end-ny*p+1:end), vik, gca);

% Row 4
nexttile;
imagesc_parts(Ef_r1, vik, gca, Up_r2, Uf_r2, Yp_r2);
ylabel('$\mathcal{Z}=[\tilde{U}_p; \tilde{U}_f; \tilde{Y}_p]_{\mathrm{opt}}$', 'Interpreter', 'latex', 'FontSize', Fig10FS);
nexttile;
imagesc_parts(Wp_r1, vik, gca, Up_r2, Uf_r2, Yp_r2);
nexttile; 
LestZopt = Yf_r1 * Ziv_opt.' / N * pinv(Wp_r1 * Ziv_opt.' / N);
imagesc_vik(LestZopt(:,1:end-ny*p), vik, gca);
nexttile;
imagesc_vik(LestZopt(:,end-ny*p+1:end), vik, gca);

Yf_iv = Yf_r1;
yr_ref = [y_ref;u_ref];

% % Transient Predictor
% tLest = zeros(ny*f,(p+f)*(nu+ny));
% for pk = p:p+f-1
%     i_k = pk-p+1;
%     
%     [~,Upk_r1,Ufk_r1] = make_Hankel(u1,pk,1);
%     [~,Ypk_r1,Yfk_r1] = make_Hankel(y1,pk,1);
%     Wpk_r1 = [Upk_r1;Ypk_r1];
%     
%     C_tKpk_hat = Yfk_r1(1:ny,:)*pinv(Wpk_r1);
%     C_tKpku_hat = C_tKpk_hat(:,1:nu*pk);
%     C_tKpky_hat = C_tKpk_hat(:,nu*pk+1:end);
%     D_hat = zeros(ny,nu);
%     tL1estk = [C_tKpku_hat D_hat zeros(ny,nu*(f-i_k)) C_tKpky_hat zeros(ny,ny*(f-i_k+1))];
%     tLest((i_k-1)*ny+1:i_k*ny,:) = tL1estk;
% end
% ImtHf = tLest(:,end-ny*f+1:end);
% tHf = eye(ny*f) - ImtHf;
% Lest = tHf\tLest(:,1:end-ny*f);
% Gf_tKpu_hat = Lest(:,1:p*nu);
% Tf_u_hat    = Lest(:,p*nu+1:(p+f)*nu);
% Gf_tKpy_hat = Lest(:,(p+f)*nu+1:(p+f)*nu+p*ny);

for kiv = 1:num_iter

% -------- get step-ahead predictor to construct approx. opt. IV ----------
% L1est = Yf_r1(1:ny,:)*Ziv.'*pinv([Up_r1;Yp_r1;Uf_r1(1:nu,:)]*Ziv.');
WZ1 = [Up_r1;Yp_r1]*Ziv.';
% L1est = [Yf_r1(1:ny,:)*Ziv.'/WZ1 zeros(ny,nu)];
% L1est = [Yf_r1(1:ny,:)*Ziv.'*WZ1.'/(WZ1*WZ1.'/N)/N zeros(ny,nu)];
L1est = [Yf_r1(1:ny,:)*Ziv.'*pinv(WZ1) zeros(ny,nu)];

CKpu_hat = L1est(:,1:p*nu);
CKpy_hat = L1est(:,p*nu+1:p*(nu+ny));
D_hat = L1est(:,end-nu+1:end);

[u_iv,y_iv] = run_pred_cl(up1,yp1,y_ref,u_ref,Cz_ss,CKpu_hat,D_hat,CKpy_hat,p,f);

% Bq = L1est(:,1:p*nu)*Mpnu+L1est(:,end-nu+1:end);
% Aq = eye(ny)-L1est(:,p*nu+1:p*(nu+ny))*Mpny;
% 
% sys_uiv = (eye(nu)+Rq\Sq/Aq*Bq)\(Rq\Tq);
% u_iv = lsim(sys_uiv,yr_ref(:,f:end).').';
% y_iv = lsim(Bq/Aq,u_iv).';

[~,Up_iv,Uf_iv] = make_Hankel(u_iv,p,fid); % p, fid
[~,Yp_iv,~] = make_Hankel(y_iv,p,fid);     % p, fid
Ziv = [Up_iv;Yp_iv;Uf_iv];
% Ziv = [Up_iv;Yp_iv];
Niv = min(size(Ziv,2),N);
Ziv = Ziv(:,1:Niv); Up_iv = Up_iv(:,1:Niv); Yp_iv = Yp_iv(:,1:Niv); Uf_iv = Uf_iv(:,1:Niv);

% estimate predictor markov parameters
WZ2 = [Up_r1;Yp_r1;Uf_r1]; WZ2 = WZ2(:,1:Niv)*Ziv.';
% Lest = Yf_r1(:,1:Niv)*Ziv.'/WZ2;
% Lest = Yf_r1(:,1:Niv)*Ziv.'*WZ2.'/(WZ2*WZ2.'/N)/N;
Lest = Yf_r1*Ziv.'*pinv(WZ2);
hat_Gf_tKp_u = get_Gf_tKp_u(Lest);
hat_Gf_tKp_y = get_Gf_tKp_y(Lest);
hat_Tf_u = get_Tf_u(Lest); % tril(get_Tf_u(Lest),-1); % enforce causality
% Lest = [hat_Gf_tKp_u hat_Gf_tKp_y hat_Tf_u];

% plotting IV data
plot(ax1_1,y_iv);
plot(ax1_2,u_iv);
nexttile(tl10); imagesc_parts(Ef_r1(:,1:Niv), vik, gca, Up_iv, Uf_iv, Yp_iv);
nexttile;       imagesc_parts(Wp_r1(:,1:Niv), vik, gca, Up_iv, Uf_iv, Yp_iv);
nexttile;       imagesc_vik([hat_Gf_tKp_u hat_Tf_u], vik, gca);
nexttile;       imagesc_vik(hat_Gf_tKp_y, vik, gca);

fro_EfZ(1+(kiv-1)*2+1) = norm(Ef_r1*Ziv.'/size(Ef_r1,2),'fro')^2;

% % ----------------------------- update IV ---------------------------------
ICuf = get_ICuf(Lest);
ICyf = get_ICyf(Lest);
ICup = get_ICup(Lest);
ICyp = get_ICyp(Lest);
ICyr = get_ICyr(Lest);
ICur = get_ICur(Lest);

Uf2Uf0 = ICuf + ICyf*hat_Tf_u;
left_Uf = eye(nu*f) - Uf2Uf0;
inv_Uf = left_Uf\[ICup+ICyf*hat_Gf_tKp_u ICyp+ICyf*hat_Gf_tKp_y ICyr ICur];
Up2Uf = inv_Uf(:,1:size(ICup,2));
Yp2Uf = inv_Uf(:,size(ICup,2)+1:size(ICup,2)+size(ICyp,2));
YRf2Uf = inv_Uf(:,end-size([ICyr ICur],2)+1:end-size(ICur,2));
URf2Uf = inv_Uf(:,end-size(ICur,2)+1:end);

% Up2Uf = get_Up2Uf(L1est);
% Yp2Uf = get_Yp2Uf(L1est);
% YRf2Uf = get_YRf2Uf(L1est);
Uf_iv = Up2Uf*Up_iv + Yp2Uf*Yp_iv + YRf2Uf*yRf + URf2Uf*uRf;        % use IV mats here?

Yf_iv0 = Yf_iv;
Yf_iv = hat_Gf_tKp_u*Up_r1 + hat_Gf_tKp_y*Yp_r1 + hat_Tf_u * Uf_iv; % use IV mats here?
Yf_iv = Yf_iv(1:ny*fid,:);
[means_yfdiff, std_devs_yfdiff] = anti_diag_stats(Yf_iv-Yf_iv0);
[means_uf, std_devs_uf] = anti_diag_stats(Uf_iv);
[means_yf, std_devs_yf] = anti_diag_stats(Yf_iv);

fro_EfUf(kiv+1) = norm(Ef_r1*Uf_iv.','fro')^2;


figure(11);
ax11_1 = subplot(2,1,1);
plot_with_shaded_std(1:length(means_uf), means_uf, std_devs_uf, [0 0.5 1]); % Blue shade
hold off;
ax11_2 = subplot(2,1,2);
plot_with_shaded_std(1:length(means_yf), means_yf, std_devs_yf, [1 0.5 0.5]); % Red shade
hold off;
linkaxes([ax11_1 ax11_2],'x');

% update IV
Ziv = R_e\[Up_iv;Yp_iv;Uf_iv];

imagesc_parts(Ef_r1,vik,fig10, R_e\Up_iv, R_e\Uf_iv, R_e\Yp_iv)
title('$E_f [U_p;U_f;Y_p]^\top$','Interpreter','latex');
fro_EfZ(kiv+1) = norm(Ef_r1*Ziv.'/size(Ef_r1,2),'fro')^2;

end

% estimate predictor markov parameters
Lest = Yf_r1*Ziv.'*pinv([Up_r1;Yp_r1;Uf_r1]*Ziv.');

% analysis - get matrix estimates
% [tGf_tKp_u_iv, tTf_u_iv, tGf_tKp_y_iv, tTf_y_iv, Gf_tKp_u_iv,Tf_u_iv,Gf_tKp_y_iv,Tf_e_iv] = L1est2mats(L1est,p,f,nu);
%%
figure();
plot(fro_EfZ); hold on
plot(fro_EfUf);
yline(fro_EfUp);
yline(fro_EfYp);
legend('Ef*Z','Ef*Uf','Ef*Up','Ef*Yp');
%% plotting
close all;
figure(1)
plot(y1.'); hold on;
plot(y_ref(1,1:Nbar+f));
plot(y_ss.');
grid on;

figure()
tl2 = tiledlayout(3,6,"TileSpacing","compact","Padding","compact");

row1_data = [Gf_tKp_u_iv Tf_u_iv Gf*tKp_u Tf_u Up2Yf_inno ];
clim_row1 = [min(row1_data,[],'all') max(row1_data,[],'all')];
nexttile(tl2,1,[1 2]);
imagesc([Gf_tKp_u_iv Tf_u_iv],clim_row1); hold on;
xline(p*nu+0.5);
ylabel("$u_p\rightarrow y_f$",'Interpreter','latex');
xlabel('est.: $\left[\Gamma_f\tilde{\mathcal{K}}_p^\mathrm{u} \;\; \mathcal{T}_f^\mathrm{u}\right]$','interpreter','latex');
nexttile(tl2,3,[1 2]);
imagesc([Gf*tKp_u Tf_u],clim_row1); hold on;
xline(p*nu+0.5);
xlabel('act.: $\left[\Gamma_f\tilde{\mathcal{K}}_p^\mathrm{u} \;\; \mathcal{T}_f^\mathrm{u}\right]$','interpreter','latex');
nexttile(tl2,5,[1 2]);
imagesc([Up2Yf_inno Tf_u],clim_row1); hold on;
xline(p*nu+0.5);
xlabel('act.: $\left[\Gamma_f(\tilde{\mathcal{K}}_p^\mathrm{u}-\tilde{A}^p\Gamma_p^\dagger\mathcal{T}_p^\mathrm{u}) \;\; \mathcal{T}_f^\mathrm{u}\right]$','Interpreter','latex');
colorbar

row2_data = [Gf_tKp_y_iv Gf*tKp_y Yp2Yf_inno];
clim_row2 = [min(row2_data,[],'all') max(row2_data,[],'all')];
nexttile(tl2,7,[1 2]);
imagesc(Gf_tKp_y_iv,clim_row2); hold on;
ylabel("$y_p\rightarrow y_f$",'Interpreter','latex');
xlabel('est.: $\Gamma_f\tilde{\mathcal{K}}_p^\mathrm{y}$','interpreter','latex');
nexttile(tl2,9,[1 2]);
imagesc(Gf*tKp_y,clim_row2); hold on;
xlabel('act.: $\Gamma_f\tilde{\mathcal{K}}_p^\mathrm{y}$','interpreter','latex');
nexttile(tl2,11,[1 2]);
imagesc(Yp2Yf_inno,clim_row2); hold on;
xlabel('act.: $\Gamma_f(\tilde{\mathcal{K}}_p^\mathrm{y}+\tilde{A}^p\Gamma_p^\dagger)$','Interpreter','latex');
colorbar

row3_data = [Tf_e_iv Tf_e];% tTf_y tTf_y_r1
clim_row3 = [min(row3_data,[],'all') max(row3_data,[],'all')];
nexttile(tl2,13,[1 3]);
imagesc(Tf_e_iv,clim_row3); hold on;
ylabel("$e_f\rightarrow y_f$",'Interpreter','latex');
xlabel('est.: $\mathcal{T}_f^\mathrm{e}$','interpreter','latex');
nexttile(tl2,16,[1 3]);
imagesc(Tf_e,clim_row3); hold on;
xlabel('act.: $\mathcal{T}_f^\mathrm{e}$','interpreter','latex');
colorbar

figure();
tl3 = tiledlayout(3,6,"TileSpacing","compact","Padding","compact");

row1_data = [tGf_tKp_u_iv tTf_u_iv tGf*tKp_u tTf_u Up2Yf_pred ];
clim_row1 = [min(row1_data,[],'all') max(row1_data,[],'all')];
nexttile(tl3,1,[1 2]);
imagesc([tGf_tKp_u_iv tTf_u_iv],clim_row1); hold on;
xline(p*nu+0.5);
ylabel("$u_p\rightarrow y_f$",'Interpreter','latex');
xlabel('est.: $\left[\tilde{\Gamma}_f\tilde{\mathcal{K}}_p^\mathrm{u} \;\; \tilde{\mathcal{T}}_f^\mathrm{u}\right]$','interpreter','latex');
nexttile(tl3,3,[1 2]);
imagesc([tGf*tKp_u tTf_u],clim_row1); hold on;
xline(p*nu+0.5);
xlabel('act.: $\left[\tilde{\Gamma}_f\tilde{\mathcal{K}}_p^\mathrm{u} \;\; \tilde{\mathcal{T}}_f^\mathrm{u}\right]$','interpreter','latex');
nexttile(tl3,5,[1 2]);
imagesc([Up2Yf_pred tTf_u],clim_row1); hold on;
xline(p*nu+0.5);
xlabel('act.: $\left[\tilde{\Gamma}_f(\tilde{\mathcal{K}}_p^\mathrm{u}-\tilde{A}^p\Gamma_p^\dagger\mathcal{T}_p^\mathrm{u}) \;\; \tilde{\mathcal{T}}_f^\mathrm{u}\right]$','Interpreter','latex');
colorbar

row2_data = [tGf_tKp_y_iv tGf*tKp_y Yp2Yf_pred];
clim_row2 = [min(row2_data,[],'all') max(row2_data,[],'all')];
nexttile(tl3,7,[1 2]);
imagesc(tGf_tKp_y_iv,clim_row2); hold on;
ylabel("$y_p\rightarrow y_f$",'Interpreter','latex');
xlabel('est.: $\tilde{\Gamma}_f\tilde{\mathcal{K}}_p^\mathrm{y}$','interpreter','latex');
nexttile(tl3,9,[1 2]);
imagesc(tGf*tKp_y,clim_row2); hold on;
xlabel('act.: $\tilde{\Gamma}_f\tilde{\mathcal{K}}_p^\mathrm{y}$','interpreter','latex');
nexttile(tl3,11,[1 2]);
imagesc(Yp2Yf_pred,clim_row2); hold on;
xlabel('act.: $\tilde{\Gamma}_f(\tilde{\mathcal{K}}_p^\mathrm{y}+\tilde{A}^p\Gamma_p^\dagger)$','Interpreter','latex');
colorbar

row3_data = [tTf_y_iv tTf_y];% tTf_y tTf_y_r1
clim_row3 = [min(row3_data,[],'all') max(row3_data,[],'all')];
nexttile(tl3,13,[1 3]);
imagesc(tTf_y_iv,clim_row3); hold on;
ylabel("$y_f\rightarrow y_f$",'Interpreter','latex');
xlabel('est.: $\tilde{\mathcal{T}}_f^\mathrm{y}$','interpreter','latex');
nexttile(tl3,16,[1 3]);
imagesc(tTf_y,clim_row3); hold on;
xlabel('act.: $\tilde{\mathcal{T}}_f^\mathrm{y}$','interpreter','latex');
colorbar

%% Helper functions

% run closed-loop simulation with real system
function [u,x,y] = run_real_cl(up,yp,x0,y_ref,u_ref,Cz,e,ny,nu,nx,Nbar,f,A,B,C,D,K)
u = nan(nu,Nbar);
y = nan(ny,Nbar);
x = nan(nx,Nbar+1);
x(:,1) = x0;

for k = 1:Nbar
    % get reference
    yrf = y_ref(:,k:k+f-1);
    urf = u_ref(:,k:k+f-1);
    
    % --------------- obtained IV ---------------
    % calculate input using initial controller
    u(:,k) = Cz(up,yp,yrf,urf);
%     u1(:,k) = rf2u*rf(:) + up2u*up1(:) + yp2u*yp1(:);
    
    % get state and output 
    x(:,k+1) = A*x(:,k) + B*u(:,k) + K*e(:,k);
    y(:,k) = C*x(:,k) + D*u(:,k) + e(:,k);

    % update up, yp
    up = [up(:,2:end) u(:,k)];
    yp = [yp(:,2:end) y(:,k)];
end
end

% run closed-loop with predictor
function [u,y] = run_pred_cl(up,yp,y_ref,u_ref,Cz_ss,CKpu,D,CKpy,p,f)
% CKpu, CKpy, D -> predictorm Markov parameter estimates
[ny,nu] = size(D);

% ---- make closed-loop system with controller & step-ahead predictor -----
[Acz,Bcz,Ccz,Dcz] = ssdata(Cz_ss);

% -> C matrix
nref = size(Cz_ss.A,1)-p*(nu+ny); % number of u & y reference states
Cpred = [CKpu CKpy zeros(ny,nref)];
Cnew = Cpred+D*Ccz;
Ccl = [Ccz;Cnew];

% -> D matrix
Dcl = [eye(nu);D]*Dcz(:,ny+1:end);

% -> B matrix
Bcl = Bcz(:,ny+1:end);

% -> A matrix
Acl = Acz;
Acl(p*(nu+ny)-ny+1:p*(nu+ny),:) = Cpred;

ss_cl = ss(Acl,Bcl,Ccl,Dcl,-1); % make system

% ---------------------- set initial conditions ---------------------------
yr_ini = y_ref(:,1:f-1);
ur_ini = u_ref(:,1:f-1);
if f > 1
    x0 = [up(:);yp(:);yr_ini(:);ur_ini(:)];
else
    x0 = [up(:);yp(:)];
end

% ------------------------ simulate system --------------------------------
uy_cl = lsim(ss_cl,[y_ref(:,f:end);u_ref(:,f:end)].',[],x0).';
u = uy_cl(1:nu,:);
y = uy_cl(nu+1:end,:);

end

% run closed-loop simulation with estimated system
function [u,y] = run_est_cl(up,yp,ref,Cz,ny,nu,Nbar,f,C_tKp_u,D,C_tKp_y)
u = nan(nu,Nbar);
y = nan(ny,Nbar);

for k = 1:Nbar
    % get reference
    rf = ref(:,k:k+f-1);
    
    % --------------- obtained IV ---------------
    % calculate input using initial controller
    u(:,k) = Cz(up,yp,rf);
%     u1(:,k) = rf2u*rf(:) + up2u*up1(:) + yp2u*yp1(:);
    
    % get output estimate
    y(:,k) = C_tKp_u*up(:) + D*u(:,k) + C_tKp_y*yp(:);

    % update up, yp
    up = [up(:,2:end) u(:,k)];
    yp = [yp(:,2:end) y(:,k)];
end

end

% make Hankel matrices with data
function [Hpf,Hp,Hf] = make_Hankel(data,s1,s2)
[ndata,Nbar] = size(data);
Hpf = mat2cell(data,ndata,ones(1,Nbar));
Hpf = Hpf(hankel(1:s1+s2,s1+s2:Nbar));
Hpf = cell2mat(Hpf);
Hp = Hpf(1:s1*ndata,:);
Hf = Hpf(s1*ndata+1:end,:);
end

% go from estimated predictor Markov parameters to full matrix estimates
function [tGf_tKpu, tTfu, tGf_tKpy, tTfy, Gf_tKpu,Tfu,Gf_tKpy,Tfe] = L1est2mats(L1est,p,f,nu)
ny = size(L1est,1);
CtKpu = L1est(:,1:p*nu);            % [C(At^p-1)Bt ... C(At)Bt C*Bt]
CtKpy = L1est(:,p*nu+1:p*(nu+ny));  % [C(At^p-1)K  ... C(At)K  C*K ]
D = L1est(:,end-nu+1:end);
Markovs_u = fliplr(mat2cell(CtKpu,ny,nu*ones(1,p))); % {C*Bt,C(At)Bt,...,C(At^p-1)Bt}
Markovs_y = fliplr(mat2cell(CtKpy,ny,ny*ones(1,p))); % {C*K ,C(At)K ,...,C(At^p-1)K }

% predictor matrix structure
toep_struct = toeplitz([p+2 ones(1,f-1)],[p+2:-1:2 ones(1,f-1)]);

% predictor form :[up,uf] -> yf
tGf_tKpu__tTfu = Markovs2BlkToeplitz(Markovs_u,D,toep_struct=toep_struct);
tGf_tKpu = tGf_tKpu__tTfu(:,1:p*nu);          % up -> yf
tTfu     = tGf_tKpu__tTfu(:,end-f*nu+1:end);  % uf -> yf

% predictor form: [yp,yf] -> yf
tGf_tKpy__tTfy = Markovs2BlkToeplitz(Markovs_y,zeros(ny,ny),toep_struct=toep_struct);
tGf_tKpy = tGf_tKpy__tTfy(:,1:end-ny*f);      % yp -> yf
tTfy     = tGf_tKpy__tTfy(:,end-ny*f+1:end);  % yf -> yf

% get matrices in innovation form
mats = (eye(ny*f)-tTfy)\[tGf_tKpu__tTfu tGf_tKpy eye(ny*f)];
Gf_tKpu = mats(:,1:p*nu);
Tfu = mats(:,p*nu+1:(p+f)*nu);
Gf_tKpy = mats(:,(p+f)*nu+1:(p+f)*nu+p*ny);
Tfe = mats(:,end-ny*f+1:end);

end

% get means and std deviations of anti-diagonals of a matrix A
function [means, std_devs] = anti_diag_stats(A)
    % Get the size of the matrix
    [m, n] = size(A);
    num_diags = m + n - 1; % Number of anti-diagonals

    % Preallocate arrays
    means = zeros(1, num_diags);
    std_devs = zeros(1, num_diags);

    % Compute statistics for each anti-diagonal
    for k = 1:num_diags
        % Extract anti-diagonal elements
        if m > 1
            anti_diag = diag(flipud(A), n - k);
        else
            anti_diag = A(k);
        end
        
        % Compute mean and standard deviation
        means(k) = mean(anti_diag);
        std_devs(k) = std(anti_diag);
    end
end

function plot_with_shaded_std(x, y, std_dev, color)
    % Ensure column vectors
    x = x(:);
    y = y(:);
    std_dev = std_dev(:);

    % Define the shaded region
    upper_bound = y + std_dev;
    lower_bound = y - std_dev;

    % Plot shaded area
    fill([x; flipud(x)], [upper_bound; flipud(lower_bound)], color, ...
        'FaceAlpha', 0.2, 'EdgeColor', 'none'); 
    hold on;

    % Plot the main line
    plot(x, y, 'Color', color, 'LineWidth', 1);

    % Formatting
    grid on;
end

% produce imagesc of correlation
function imagesc_parts(Ef,vik,axh,varargin)
narginchk(4,inf);
Z = [];
nvarargin = length(varargin);
xvals = zeros(1,nvarargin);
for k = 1:nvarargin
    Z = [Z;varargin{k}];
    xvals(k) = size(varargin{k},1);
end
corEfZ = Ef*Z.'/size(Ef,2);
imagesc_vik(corEfZ,vik,axh);
xvals = cumsum(xvals)+0.5;
for k = 1:nvarargin-1
    xline(xvals(k),'LineWidth',1.5);
end
end

% produce colour-scaled imagesc using vik
function imagesc_vik(CorMat,vik,axh)
imagesc(axh,CorMat);
c_absmax = max(abs(CorMat),[],'all');
colormap(vik);
colorbar;
caxis([-c_absmax c_absmax]);
end