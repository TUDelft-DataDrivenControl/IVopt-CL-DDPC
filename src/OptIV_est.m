%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           DDPC using an Optimal-IV
%           Authors: R. Dinkla, T. Oomen, J.W. van Wingerden
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solves a constrained least-squares regression problem
clear;
close all;
rng default;

addpath('..\bin');

%% problem parameters
f = 10;
fid = 1;
p = 20;
ny = 1;
nu = 1;
N = 1e4;
Nbar = p + fid + N -1;

% controller weights
dRk= 1*eye(nu);  dR = kron(speye(f),dRk);
Rk = 0*eye(nu);   R = kron(speye(f),Rk);
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

% --------------------- create actual matrices ----------------------------
At = A-K*C;
Bt = B-K*D;

pf_max = max(p,f);

% Yf =  Gf * X0 +  Tf_u * Uf +              Tf_e * Ef
% Yf = tGf * X0 + tTf_u * Uf + tTf_y * Yf +        Ef

% make extended observability matrices
G = make_ext_obsv(C,A,pf_max);
Gp = G(1:p*ny,:);
Gf = G(1:f*ny,:); clear G;
tGf = make_ext_obsv(C,At,f);

% make extended reversed controllability matrices
tKp_u = make_ext_rev_ctrb(At,Bt,p);
tKp_y = make_ext_rev_ctrb(At,K,p);

% make block-lower toeplitz matrices
T_u = make_blk_tril_toeplitz(A,B,C,D,pf_max);
Tf_u = T_u(1:f*ny,1:f*nu); 
Tp_u = T_u(1:p*ny,1:p*nu); clear T_u;
tT_u= make_blk_tril_toeplitz(At,Bt,C,D,pf_max);
tTf_u = tT_u(1:f*ny,1:f*nu); 
tTp_u = tT_u(1:p*ny,1:p*nu); clear tT_u;
tTf_y = make_blk_tril_toeplitz(At,K,C,zeros(ny,ny),f);
T_e = make_blk_tril_toeplitz(A,K,C,eye(ny,ny),pf_max);
Tp_e = T_e(1:p*ny,1:p*ny);
Tf_e = T_e(1:f*ny,1:f*ny); clear T_e;
% note: inv(I - tTf_y) * [tGf tTf_u eye(ny*f)] == [Gf Tf_u Tf_e]
% (eye(ny*f)-tTf_y)\[tGf tTf_u eye(ny*f)]-[Gf Tf_u Tf_e] % <- check ==0

% useful constants
Atp_pinvGp = At^p*pinv(Gp);
Up2X0 = tKp_u-Atp_pinvGp*Tp_u; % Up -> X0
Yp2X0 = tKp_y+Atp_pinvGp;      % Yp -> X0
Ep2X0 = -Atp_pinvGp*Tp_e;      % Ep -> X0

% predictor form matrices:
Up2Yf_pred = tGf*Up2X0;
Uf2Yf_pred = tTf_u;
Yp2Yf_pred = tGf*Yp2X0;
Yf2Yf_pred = tTf_y;
Ep2Yf_pred = tGf*Ep2X0;
Ef2Yf_pred = eye(ny*f,ny*f);

% innovation form matrices:
Up2Yf_inno = Gf*Up2X0;
Uf2Yf_inno = Tf_u;
Yp2Yf_inno = Gf*Yp2X0;
Ep2Yf_inno = Gf*Ep2X0;
Ef2Yf_inno = Tf_e;

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
[u1,x1,y1] = run_real_cl(up1,yp1,x1_0,ref,Cz0,e1,ny,nu,nx,Nbar,f,A,B,C,D,K);

% creating Hankel matrices
[~,Up_r1,Uf_r1] = make_Hankel(u1,p,fid);
[~,Yp_r1,Yf_r1] = make_Hankel(y1,p,fid);
[~,Ep_r1,Ef_r1] = make_Hankel(e1,p,fid);
[~,Rp,Rf] = make_Hankel(ref(:,1:Nbar),p,fid);

%% -------------------- Get system estimate using IV -----------------------
% creating IV matrix
Ziv = [Up_r1;Yp_r1;Rf];

for kiv = 1:5
figure(10);
imagesc(Ef_r1*Ziv.');
colorbar;

% estimate predictor markov parameters
L1est = Yf_r1*Ziv.'*pinv([Up_r1;Yp_r1;Uf_r1]*Ziv.');
CtKp_u_iv = L1est(:,1:p*nu);
D_iv = L1est(:,p*nu+1:(p+1)*nu);
CtKp_y_iv = L1est(:,end-p*ny+1:end);

% analysis - get matrix estimates
% [tGf_tKp_u_iv, tTf_u_iv, tGf_tKp_y_iv, tTf_y_iv, Gf_tKp_u_iv,Tf_u_iv,Gf_tKp_y_iv,Tf_e_iv] = L1est2mats(L1est,p,f,nu);

% analysis - error
%

% bootstrapping: iterative optimal IV estimation
[u_iv,y_iv] = run_est_cl(up1,yp1,ref,Cz0,ny,nu,Nbar,f,CtKp_u_iv,D_iv,CtKp_y_iv);

figure(2);
plot(y_iv);

% creating Hankel matrices
[~,Up_iv,Uf_iv] = make_Hankel(u_iv,p,fid);
[~,Yp_iv,Yf_iv] = make_Hankel(y_iv,p,fid);

% update IV
Ziv = [Up_iv;Yp_iv;Uf_iv];
end

% estimate predictor markov parameters
L1est = Yf_r1*Ziv.'*pinv([Up_r1;Yp_r1;Uf_r1]*Ziv.');

% analysis - get matrix estimates
[tGf_tKp_u_iv, tTf_u_iv, tGf_tKp_y_iv, tTf_y_iv, Gf_tKp_u_iv,Tf_u_iv,Gf_tKp_y_iv,Tf_e_iv] = L1est2mats(L1est,p,f,nu);

%% plotting
close all;
figure(1)
plot(y1.'); hold on;
plot(ref(1,1:Nbar+f));
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
function [u,x,y] = run_real_cl(up,yp,x0,ref,Cz,e,ny,nu,nx,Nbar,f,A,B,C,D,K)
u = nan(nu,Nbar);
y = nan(ny,Nbar);
x = nan(nx,Nbar+1);
x(:,1) = x0;

for k = 1:Nbar
    % get reference
    rf = ref(:,k:k+f-1);
    
    % --------------- obtained IV ---------------
    % calculate input using initial controller
    u(:,k) = Cz(up,yp,rf);
%     u1(:,k) = rf2u*rf(:) + up2u*up1(:) + yp2u*yp1(:);
    
    % get state and output 
    x(:,k+1) = A*x(:,k) + B*u(:,k) + K*e(:,k);
    y(:,k) = C*x(:,k) + D*u(:,k) + e(:,k);

    % update up, yp
    up = [up(:,2:end) u(:,k)];
    yp = [yp(:,2:end) y(:,k)];
end
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

% make extended observability matrix
function Gs = make_ext_obsv(C,A,s)
Gs = cell(s,1);
Gs{1} = C;
for k = 2:s
    Gs{k}  = Gs{k-1}*A;
end
Gs  = cell2mat(Gs);
end

% make extended reversed controllability matrix
function Ks = make_ext_rev_ctrb(A,B,s)
Ks = cell(1,s);
Ks{s} = B;
for k = s:-1:2
     Ks{k-1} = A  *  Ks{k};
end
Ks = cell2mat(Ks);
end

% make block-lower toeplitz matrix from Markov parameters
function BlkTrilToep = Markovs2BlkToeplitz(Markovs,D,opts)
arguments
    Markovs cell
    D       double
    opts.depth (1,1) double
    opts.toep_struct = [];
end
Markovs2 = [{zeros(size(D))}, {D}, Markovs(:).'];
if isfield(opts,'depth')
    s = opts.depth;
    toep_struct = toeplitz(2:s+1,[2,ones(1,s-1)]);
elseif ~isempty(opts.toep_struct)
    toep_struct = opts.toep_struct;
else
    error('Provide either number of block-rows or toeplitz structure');
end
BlkTrilToep = cell2mat(Markovs2(toep_struct));
end

% make block-lower toeplitz matrix
function BlkTrilToep = make_blk_tril_toeplitz(A,B,C,D,s)
Markovs = get_Markovs(A,B,C,s-1);
BlkTrilToep = Markovs2BlkToeplitz(Markovs,D,depth=s);
end

% get Markov parameters: {CB, CAB, C(A^2)B, ..., C(A^(s-1))B}
function Markovs = get_Markovs(A,B,C,s)
Markovs = cell(1,s);
Markovs{1} = B;
for k = 2:s
     Markovs{k} = A  *  Markovs{k-1};
end
Markovs = cellfun(@(x) C*x, Markovs, UniformOutput=false);
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