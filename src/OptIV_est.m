%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           DDPC using an Optimal-IV
%           Authors: R. Dinkla, T. Oomen, J.W. van Wingerden
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solves a constrained least-squares regression problem
clear;
close all;
rng default;

addpath(genpath('..\bin'));
addpath(genpath('.\fun'));
% addpath(genpath("C:\Users\rogierdinkla\Documents\MATLAB"));
% addpath(genpath("C:\Users\rogierdinkla\Documents\CasADi\casadi-3.6.7-windows64-matlab2018b"));
% % color maps: https://www.fabiocrameri.ch/colourmaps/
% load('vik.mat');

%% problem parameters
f = 15;
fid = f;
p = 20;
ny = 1;
nu = 1;
N = 1e4;%100*(p*(nu+ny)+f*nu);
Nbar = p + fid + N -1;

% controller weights
dRk= 1*eye(nu);  dR = kron(speye(f),dRk);
Rk = 1*eye(nu);   R = kron(speye(f),Rk);
Qk = 1e2*eye(ny); Q = kron(speye(f),Qk);

% reference
Nbar_ref = ceil((Nbar-fid)/f)*f+2*f-1;
y_ref = idinput(Nbar_ref,'prbs',[],[-1 1]).';
% y_ref = sign(sin((1:3:3*(Nbar-fid+2*f))*2*pi/(Nbar/20)));
% y_ref = 0.5*sign(sin(pi/2+(1:3:3*(Nbar-fid+2*f))*2*pi/(Nbar/60))) + ...
%         0.5*sign(sin((1:3:3*(Nbar-fid+2*f))*2*pi/(Nbar/60)));
% y_ref = y_ref(:,1:end-1);
y_ref = repmat(y_ref,ny,1);

% get solver: uf_ref, yf_ref, up, yp -> u_k
[up2usol, yp2usol, urf2usol, yrf2usol, get_ICuf, get_ICyf, get_ICup,...
    get_ICyp, get_ICyr, get_ICur, calc_u_v2] = get_solver(nu,ny,p,f,Q,R,dR);

%% ============================= make system ==============================
% make nominal system
% nx = 2;
% sys = drss(nx,ny,nu);% sys = c2d(sys,1e-1);
% [A,B,C,D] = ssdata(sys); D = 0*D; sys = ss(A,B,C,D,-1);
% model_Favoreel1999; % D = 0...
% model_DoubleTank; R_e = 1e0*Re;
[plant,nx,nu,ny,A,B,C,D,K,R_e] = model_Landau1995(Re=4.81);
% [plant,Cz,nx,nu,ny,A,B,C,D,K,Re] = model_Wang2023(Re=5);
% R_e = Re;
% [~,Vc,Uc] = lncf(Cz); % Cz = Vc\Uc

% set noise properties
% R_e = 1e-0*tril(rand(ny)*2);
% if any(diag(R_e)==0) % ensure Re > 0
%     var_e_diag = diag(R_e);
%     var_e_diag(var_e_diag==0) = 1e-2;
%     indiag = 1:ny;
%     R_e(indiag + (indiag - 1)*ny) = var_e_diag;
% end
% R_e = R_e*R_e.';
% K = place(A.',C.',linspace(0.9,0.95,nx)).';

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
ySF = diag(mean(abs(y1),2));
u = uSF\u; u1 = uSF\u1; u2 = uSF\u2; up1 = uSF\up1; u_ref = uSF\u_ref;
y = ySF\y; y1 = ySF\y1; y2 = ySF\y2; yp1 = ySF\yp1; y_ref = ySF\y_ref;
e = ySF\e; e1 = ySF\e1; R_e = ySF.^2\R_e;
B = B*uSF; D = ySF\D*uSF; C = ySF\C; K = K*ySF;
sys2 = ss(A,[B K],C,[D eye(ny)],-1);

% creating Hankel matrices
[~,Up_r1,Uf_r1] = make_Hankel(u1,p,fid); [~,Up_r2,Uf_r2] = make_Hankel(u2,p,fid);
[~,Yp_r1,Yf_r1] = make_Hankel(y1,p,fid); [~,Yp_r2,Yf_r2] = make_Hankel(y2,p,fid); 
[~,Ep_r1,Ef_r1] = make_Hankel(e1,p,fid);
[~,yuRp,yuRf] = make_Hankel([y_ref;u_ref],p,2*f-1); yuRp = yuRp(:,1:N); yuRf(:,1:N);
[~,yRp,yRf] = make_Hankel(y_ref,p,2*f-1); yRp = yRp(:,1:N); yRf = yRf(:,1:N);
[~,uRp,uRf] = make_Hankel(u_ref,p,2*f-1); uRp = uRp(:,1:N); uRf = uRf(:,1:N);
Wp_r1 = [Up_r1;Uf_r1;Yp_r1];

%% make actual matries and scaled controller components

% build matrix representation of closed-loop
up2uk  = uSF\up2usol(Lest_r0) *kron(speye(p),uSF);
yp2uk  = uSF\yp2usol(Lest_r0) *kron(speye(p),ySF);
yrf2uk = uSF\yrf2usol(Lest_r0)*kron(speye(f),ySF);
urf2uk = uSF\urf2usol(Lest_r0)*kron(speye(f),uSF);

ICuf = kron(speye(f),uSF)\get_ICuf(Lest_r0)*kron(speye(f),uSF);
ICyf = kron(speye(f),uSF)\get_ICyf(Lest_r0)*kron(speye(f),ySF);
ICup = kron(speye(f),uSF)\get_ICup(Lest_r0)*kron(speye(p),uSF);
ICyp = kron(speye(f),uSF)\get_ICyp(Lest_r0)*kron(speye(p),ySF);
ICyr = kron(speye(f),uSF)\get_ICyr(Lest_r0)*kron(speye(2*f-1),ySF);
ICur = kron(speye(f),uSF)\get_ICur(Lest_r0)*kron(speye(2*f-1),uSF);

% --------------------- create actual matrices ----------------------------
[Up2Yf_inno, Yp2Yf_inno, Uf2Yf_inno] = get_actual_matrices(A,B,C,D,K,p,f);

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

ax1_1.LineStyleOrderIndex = ax1_1.LineStyleOrderIndex; % [1]
ax1_1.LineStyleOrder = {'-o','-+','-*','-x','-s','-d','-v','->','-h','-^'};
ax1_1.ColorOrder = [1 0 0; 0 1 0; 0 0 1; 0 1 1; 1 0 1];
hold(ax1_1,'on') % [2]
ax1_2.LineStyleOrderIndex = ax1_2.LineStyleOrderIndex; % [1]
ax1_2.LineStyleOrder = ax1_1.LineStyleOrder;
ax1_2.ColorOrder = ax1_1.ColorOrder;
hold(ax1_2,'on') % [2]
% [1] This line forces the axes to behave consistently between releases before and after r2019a
% [2] Important to preserve the property orders that were just set.

% -------------------- Get system estimate using IV -----------------------
num_iter = 3;

% for plotting RMSE to optimal IV
fig3 = figure(3);
RMSE_Yiv = nan(1,num_iter);
RMSE_Uiv = RMSE_Yiv;

% optimal IV matrix
Ziv_opt = [Up_r2;Uf_r2;Yp_r2];

% creating initial IV matrix
Ziv = [yRp;yRf;Yp_r1]; %[Up_r2;Yp_r2;Uf_r2];
% Ziv = [Up_r1;Uf_r1;Yp_r1];

% plotting correlations and estimate of system parameters
plot_corrs;

% LHS = (eye(nu*f)-ICuf)*Uf_r1-ICyr*yRf-ICur*uRf-ICup*Up_r1-ICyp*Yp_r1;
% LHS = pinv([ICyf;eye(ny*f)])*[LHS;Yf_r1];

for kiv = 1:num_iter
% if kiv < 3

% using multiple-step with controller information
%
% ================================ CL-SPC =================================
% WZ1 = [Up_r1;Yp_r1]*Ziv.';
% L1est = [Yf_r1(1:ny,:)*Ziv.'*pinv(WZ1) zeros(ny,nu)];
% % L1est = [LHS(1:ny,:)*Ziv.'*pinv(WZ1) zeros(ny,nu)];
% C_tKpu_hat = L1est(:,1:p*nu);
% C_tKpy_hat = L1est(:,p*nu+1:p*(nu+ny));
% D_hat = L1est(:,end-nu+1:end);
% 
% % construct multiple-step-ahead predictor
% tLest_u = zeros(ny*f,nu*(p+f));
% tLest_u(1:ny,1:nu*(p+1)) = [C_tKpu_hat D_hat];
% tLest_y = zeros(ny*f,ny*(p+f));
% tLest_y(1:ny,1:ny*p) = C_tKpy_hat;
% for k = 2:f
%     tLest_u((k-1)*ny+1:k*ny,:) = circshift(tLest_u((k-2)*ny+1:(k-1)*ny,:),nu,2);
%     tLest_y((k-1)*ny+1:k*ny,:) = circshift(tLest_y((k-2)*ny+1:(k-1)*ny,:),ny,2);
% end
% tHf = eye(ny*f)-tLest_y(:,end-ny*f+1:end);
% Lest = tHf\[tLest_u(:,1:p*nu) tLest_y(:,1:p*ny) tLest_u(:,end-nu*f+1:end)];
% hat_Gf_tKp_u = Lest(:,1:p*nu);
% hat_Gf_tKp_y = Lest(:,p*nu+1:p*(nu+ny));
% hat_Tf_u = Lest(:,end-nu*f+1:end);
%{%
% ================================ SPC ====================================
WpZ = [Up_r1;Yp_r1;Uf_r1]*Ziv.';
Lest = Yf_r1*Ziv.'*pinv(WpZ);
% Lest = LHS*Ziv.'*pinv(WpZ);
cond(WpZ)
hat_Gf_tKp_u = Lest(:,1:p*nu);
hat_Gf_tKp_y = Lest(:,p*nu+1:p*(nu+ny));
hat_Tf_u = Lest(:,end-nu*f+1:end);
% hat_Tf_u = tril(hat_Tf_u,-1);

HfEf = Yf_r1 - Lest*[Up_r1;Yp_r1;Uf_r1];
[Ef_hat,hat_Hf] = qr(HfEf.','econ'); hat_Hf = hat_Hf.'; Ef_hat = Ef_hat.';

% Step 2: Normalize diagonal blocks of R to identity
D_block = cell(f,1);
for i = 1:f
    row_idx = (i-1)*ny + 1:i*ny;
    D_block{i} = hat_Hf(row_idx, row_idx);  % Diagonal ny-by-ny block
end
D_blocks = sparse(blkdiag(D_block{:}));
hat_Hf = hat_Hf/D_blocks;
Ef_hat = D_blocks*Ef_hat;
tHf = inv(hat_Hf);

% ========================= Transient Predictor ===========================
% tLest = zeros(ny*f,(p+f)*(nu+ny));
% for pk = p:p+f-1
%     i_k = pk-p+1;
% 
%     [~,Upk_r1,Ufk_r1] = make_Hankel(u1,pk,1); Upk_r1 = Upk_r1(:,1:N); Ufk_r1 = Ufk_r1(:,1:N);
%     [~,Ypk_r1,Yfk_r1] = make_Hankel(y1,pk,1); Ypk_r1 = Ypk_r1(:,1:N); Yfk_r1 = Yfk_r1(:,1:N);
%     Wpk_r1 = [Upk_r1;Ypk_r1]*Ziv.';
%     C_tKpk_hat = Yfk_r1(1:ny,:)*Ziv.'*pinv(Wpk_r1);
%     C_tKpku_hat = C_tKpk_hat(:,1:nu*pk);
%     C_tKpky_hat = C_tKpk_hat(:,nu*pk+1:end);
%     D_hat = zeros(ny,nu);
%     tL1estk = [C_tKpku_hat D_hat zeros(ny,nu*(f-i_k)) C_tKpky_hat zeros(ny,ny*(f-i_k+1))];
%     tLest((i_k-1)*ny+1:i_k*ny,:) = tL1estk;
% end
% ImtHf = tLest(:,end-ny*f+1:end);
% tHf = eye(ny*f) - ImtHf;
% Lest = tHf\tLest(:,1:end-ny*f);
% hat_Gf_tKpu = Lest(:,1:p*nu);
% hat_Tf_u    = Lest(:,p*nu+1:(p+f)*nu);
% hat_Gf_tKpy = Lest(:,(p+f)*nu+1:(p+f)*nu+p*ny);

% end

% build matrix representation of closed-loop
yrf_x0 = y_ref(:,1:2*f-1); yrf = y_ref(:,2*f:end); yrf = reshape(yrf,ny*f,[]);
urf_x0 = u_ref(:,1:2*f-1); urf = u_ref(:,2*f:end); urf = reshape(urf,nu*f,[]);
[sys_cl,x0_cl] = get_pred_fstep_cl(f,p,ny,nu,ICuf,ICyf,ICup,ICyp,ICyr,ICur,hat_Tf_u,hat_Gf_tKp_u,hat_Gf_tKp_y,up1,yp1,yrf_x0,urf_x0);
x0_cl = 0*x0_cl;
% cl_data = iddata([reshape(u1(:,1:f*size(urf,2)),nu*f,[]);reshape(y1(:,1:f*size(yrf,2)),ny*f,[])].',[urf;yrf].');
% x0_cl = findstates(idss(sys_cl),cl_data);
[uy_iv,~,x_iv] = lsim(sys_cl,[urf;yrf].',[],x0_cl); uy_iv = uy_iv.'; x_iv = x_iv.';
x_iv_new = sys_cl.A*x_iv(:,end) + sys_cl.B*[urf(:,end);yrf(:,end)];
uy_iv = [uy_iv sys_cl.C*x_iv_new];
u_iv = reshape(uy_iv(1:f*nu,:),nu,[]);     u_iv = u_iv(:,1:Nbar);
y_iv = reshape(uy_iv(f*nu+1:end,:),ny,[]); y_iv = y_iv(:,1:Nbar);

% analysis A matrix of closed-loop system
[A,B,C,D] = ssdata(sys_cl);
plotABCD;
%}

% ----------------- identifying closed-loop for IV ------------------------
%{
% if kiv == 1
    L_cl = [Uf_r1;Yf_r1]*pinv([yRp;yRf]); %TODO: change... uRp, uRf make sense...? see line 344
% else
    % L_cl = 1/size(yRp,2)*[Uf_r1;Yf_r1]*Ziv.'*pinv(1/size(yRp,2)*[yRp;yRf]*Ziv.');%[Up_iv;Uf_iv;Yf_iv]
% end
L_clu = L_cl(1:nu*f,:);
L_cly = L_cl(nu*f+1:end,:);
L_clu = L_clu.*kron(tril(ones(f,p+2*f-1),p+f-1),ones(ny,nu));
L_cly = L_cly.*kron(tril(ones(f,p+2*f-1),p+f-1),ones(ny,ny));
L_clu = blk_toeplitz_mean(L_clu,ny,nu);
L_cly = blk_toeplitz_mean(L_cly,ny,ny);
L_cl = [L_clu; L_cly];

y_ref2 = [zeros(ny,p) y_ref];
[~,yRp2,yRf2] = make_Hankel(y_ref2,p,2*f-1);
UY_iv = L_cl*[yRp2;yRf2];
U_iv = UY_iv(1:nu*f,:);
U_iv = flipud(blk_toeplitz_mean(flipud(U_iv),nu,1));
% u_iv = diag_stats(flipud(U_iv));
u_iv = U_iv(:,1); u_iv = [reshape(u_iv,nu,[]),U_iv(end-nu+1:end,2:end)];
Y_iv = UY_iv(nu*f+1:end,:);
Y_iv = flipud(blk_toeplitz_mean(flipud(Y_iv),ny,1));
% y_iv = diag_stats(flipud(Y_iv));
y_iv = Y_iv(:,1); y_iv = [reshape(y_iv,ny,[]),Y_iv(end-ny+1:end,2:end)];
u_iv = u_iv(:,1:length(u1));
y_iv = y_iv(:,1:length(y1));
%}
% -------------------------------------------------------------------------

% analysis - rmse w.r.t optimal IVs
RMSE_Yiv(kiv) = rms(y_iv-y2(1:length(y_iv)));
RMSE_Uiv(kiv) = rms(u_iv-u2(1:length(u_iv)));
figure(fig3);
subplot(2,1,1); grid on;
plot(RMSE_Yiv,'-+'); ylabel('RMSE_{yiv}')
subplot(2,1,2); grid on;
plot(RMSE_Uiv,'-+'); ylabel('RMSE_{uiv}')
xlabel('iteration');

% TODO: having reflected unstable eigenvalues of closed-loop -> update
% prediction

% update IV
[~,Up_iv,Uf_iv] = make_Hankel(u_iv,p,fid);
[~,Yp_iv,Yf_iv] = make_Hankel(y_iv,p,fid);
Ziv = [Up_iv;Yp_iv;Uf_iv];

% new estimate
WZ1 = [Up_r1;Uf_r1;Yp_r1];%*Ziv.';
% tLest = tHf*Yf_r1 * Ziv.' / N * pinv(WZ1 * Ziv.' / N);
% Lest = tHf\tLest;
Lest = Yf_r1 * Ziv.' / N * pinv(WZ1 * Ziv.' / N);
hat_Gf_tKp_u = Lest(:,1:p*nu);
hat_Gf_tKp_y = Lest(:,end-p*ny+1:end);
hat_Tf_u = Lest(:,p*nu+1:(p+f)*nu);

% plotting IV data
plot(ax1_1,y_iv); ylim(ax1_1,[-2 2]);
plot(ax1_2,u_iv); ylim(ax1_2,[-10 10]);
nexttile(tl10); imagesc_parts(Ef_r1, Up_iv, Uf_iv, Yp_iv);
nexttile;       imagesc_parts(Wp_r1, Up_iv, Uf_iv, Yp_iv);
nexttile;       imagesc_vik([hat_Gf_tKp_u hat_Tf_u]);
nexttile;       imagesc_vik(hat_Gf_tKp_y);

hat_Tf_u = tril(hat_Tf_u,-1);

end

% estimate predictor markov parameters
Lest = Yf_r1*Ziv.'*pinv([Up_r1;Yp_r1;Uf_r1]*Ziv.');

% analysis - get matrix estimates
% [tGf_tKp_u_iv, tTf_u_iv, tGf_tKp_y_iv, tTf_y_iv, Gf_tKp_u_iv,Tf_u_iv,Gf_tKp_y_iv,Tf_e_iv] = L1est2mats(L1est,p,f,nu);

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

% run closed-loop with single-step ahead predictor
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
% function [Hpf,Hp,Hf] = make_Hankel(data,s1,s2)
% [ndata,Nbar] = size(data);
% Hpf = mat2cell(data,ndata,ones(1,Nbar));
% Hpf = Hpf(hankel(1:s1+s2,s1+s2:Nbar));
% Hpf = cell2mat(Hpf);
% Hp = Hpf(1:s1*ndata,:);
% Hf = Hpf(s1*ndata+1:end,:);
% end
function [Hpf, Hp, Hf] = make_Hankel(data, s1, s2)
% Efficiently constructs past-future Hankel matrices from input data
% Inputs:
%   data: (ndata x Nbar) time series data matrix
%   s1: number of past block rows
%   s2: number of future block rows
% Outputs:
%   Hpf: (ndata*(s1+s2) x (Nbar - s1 - s2 + 1)) full Hankel matrix
%   Hp:  (ndata*s1 x num_cols) past block
%   Hf:  (ndata*s2 x num_cols) future block

[ndata, Nbar] = size(data);
num_cols = Nbar - s1 - s2 + 1;

Hpf = zeros(ndata * (s1 + s2), num_cols);

for i = 1:(s1 + s2)
    Hpf((i-1)*ndata+1:i*ndata, :) = data(:, i:i+num_cols-1);
end

Hp = Hpf(1:ndata*s1, :);
Hf = Hpf(ndata*s1+1:end, :);
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
function imagesc_parts(Ef,varargin)
narginchk(2,inf);
Z = [];
nvarargin = length(varargin);
xvals = zeros(1,nvarargin);
for k = 1:nvarargin
    Z = [Z;varargin{k}];
    xvals(k) = size(varargin{k},1);
end
corEfZ = Ef*Z.'/size(Ef,2);
imagesc_vik(corEfZ);
xvals = cumsum(xvals)+0.5;
for k = 1:nvarargin-1
    xline(xvals(k),'LineWidth',1.5);
end
end

function [ad_mean, ad_var] = hankel_anti_diag_stats(H)
% Computes the mean and variance of the anti-diagonals of a Hankel matrix H

    [m, n] = size(H);
    num_antis = m + n - 1;
    
    ad_mean = zeros(1, num_antis);
    ad_var = zeros(1, num_antis);
    
    for k = 1:num_antis
        % Get indices (i, j) such that i + j == k (constant anti-diagonals)
        vals = H((1:m)'+(0:n-1) == k);
        
        ad_mean(k) = mean(vals);
        ad_var(k) = var(vals);
    end
end