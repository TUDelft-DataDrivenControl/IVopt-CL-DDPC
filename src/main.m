%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           DDPC using an Optimal-IV
%           Authors: R. Dinkla, T. Oomen, J.W. van Wingerden
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
close all;
clc;
rng default;

addpath('..\bin');
addpath(genpath("C:\Users\rogierdinkla\Documents\MATLAB"));
addpath(genpath("C:\Users\rogierdinkla\Documents\CasADi\casadi-3.6.7-windows64-matlab2018b"));
% color maps: https://www.fabiocrameri.ch/colourmaps/
load('vik.mat');

%% Simulation settings
% plant model
% model_DoubleTank; plant.Ts = 1; R_e = 1e-2;
% model_Bemporad2002; plant.Ts = 1; K = acker(A.',C.',zeros(nx,1)).'; R_e = 1e2;
[plant,nx,nu,ny,A,B,C,D,K,R_e] = model_Landau1995(Re=4.81e-1); %K = acker(A.',C.',zeros(nx,1)).';
plant.u = {'u','e'};
plant.y = 'y';

% ----------------------- make initial controller -------------------------
W1 = makeweight(33,5,0.5);  W1 = c2d(W1,plant.Ts,'tustin');
W3 = makeweight(0.5,20,20); W3 = c2d(W3,plant.Ts,'tustin');%plant(:,2);
[Cz,CL,gamma,info] = mixsyn(plant(:,1:nu),W1,[],W3);
Cz.u = 'error';
Cz.y = 'u';

% L = plant(:,1:nu)*Cz;
% I = eye(size(L));
% S = feedback(I,L); 
% T= I-S;
% 
% figure();
% sigma(S,'b',1/W1,'b--',T,'r',1/W3,'r--',{0.1,1000})
% legend('S','1/W1','T','1/W3');

% ---------------------- set SPC controller settings ----------------------
% window lengths
f = 20;
fid = f;
p = 20;
N = 1e5;

% weights
dRk= 1*eye(nu);  dR = kron(speye(f),dRk);
Rk = 1*eye(nu);   R = kron(speye(f),Rk);
Qk = 1e2*eye(ny); Q = kron(speye(f),Qk);

% ----------------------- set simulation lengths --------------------------
Nbar = p + fid + N -1; % sim. length of initial controller
Ncl = 1500;            % sim. length of SPC

% ------------------------------ set references ---------------------------
Ncl_ref = ceil((Ncl-fid)/f)*f+2*f-1;
y_ref = make_reference(Nbar+Ncl_ref,ny);%,var_dyr=0.01);
% y_ref = idinput(Nbar+Ncl_ref,'prbs',[],[-1 1]).';
% y_ref = repmat(y_ref,ny,1);
[~,Rp_r0,Rf_r0] = make_Hankel(y_ref(:,1:Nbar),p,f);

% -> references for SPC controllers
yr = y_ref(:,end-Ncl_ref+1:end);        % y-ref
P0 = dcgain(plant(:,1:nu));             % DC gain
ur = P0\yr;                             % u-ref

% ----------------------- make initial SPC controller -------------------------
% get_solver;
% e0 = mvnrnd(zeros(ny,1),R_e,Nbar).'; % create innovation noise
% u0 = mvnrnd(zeros(nu,1),1,Nbar).'; % create open-loop inputs
% [y0,~,x0] = lsim(plant,[u0;e0].',[]); y0 = y0.'; x0 = x0.';
% [~,Up_r0,Uf_r0] = make_Hankel(u1,p,f); % - data w/ noise
% [~,Yp_r0,Yf_r0] = make_Hankel(y1,p,f);
% Lf0 = Yf_r0*pinv([Up_r0;Yp_r0;Uf_r0]);
% Cz = Lf_2_SPC(Lf0,up2usol,yp2usol,urf2usol,yrf2usol,nu,ny,p,f);


%% Data generation (closed-loop)

% make closed-loop system
fbsum = sumblk('error = r - y');
Tcl = connect(Cz,plant,fbsum,{'r','e'},{'u','y'});

% ------------------------ simulation with noise --------------------------
e = mvnrnd(zeros(ny,1),R_e,Nbar+Ncl).'; % create innovation noise

% simulate system with noise
[uy1,~,x1] = lsim(Tcl,[y_ref(:,1:Nbar);e(:,1:Nbar)],[]); uy1 = uy1.';
u1 = uy1(1:nu,:); y1 = uy1(nu+1:end,:); clear uy1;       % get inputs and outputs
[~,~,xc1] = lsim(Cz,y_ref(:,1:Nbar)-y1,[]); xc1 = xc1.'; % get controller states

% create Hankel matrices
[~,Up_r1,Uf_r1] = make_Hankel(u1,p,f); % - data w/ noise
[~,Yp_r1,Yf_r1] = make_Hankel(y1,p,f);

% ---------------------- simulation without noise -------------------------
[uy2,~,x2] = lsim(Tcl,[y_ref(:,1:Nbar);zeros(ny,Nbar)],[]); uy2 = uy2.';
u2 = uy2(1:nu,:); y2 = uy2(nu+1:end,:); clear uy2;

% create Hankel matrices
[~,Up_r2,Uf_r2] = make_Hankel(u2,p,f); % - data w/o noise
[~,Yp_r2,Yf_r2] = make_Hankel(y2,p,f); 

%% Controller settings
% get_solver;

% for initial state of SPC controllers
up1 = u1(:,end-p+1:end); up1 = up1(:);  % past input data
yp1 = y1(:,end-p+1:end); yp1 = yp1(:);  % past output data
urf = ur(:,1:f); urf = urf(:);          % future input references
yrf = yr(:,1:f); yrf = yrf(:);          % future output references

% ----------------------------- get IVs -----------------------------------
% 1) optimal IV
Zopt = [Up_r2;Yp_r2;Uf_r2];

% 2) open-loop IV (i.e. least-squares regression)
W1 = [Up_r1;Yp_r1;Uf_r1];

% 3) composed IV using LCF
[frac,Vc,Uc] = lncf(ss(Cz));
Hcv = make_blk_tril_toeplitz(Vc.A,Vc.B,Vc.C,Vc.D,f);
Hcu = make_blk_tril_toeplitz(Uc.A,Uc.B,Uc.C,Uc.D,f);
IV_Theta = Hcv*Uf_r1 + Hcu*Yf_r1;
Ziv3 = [Up_r1; Yp_r1; IV_Theta; Rf_r0];

% 4) approximation of optimal IV w/o controller information
[u_iv4,y_iv4] = approx_IV_no_controller_info(u1,y1,y_ref(:,1:Nbar),p,p,f);
[~,Up_iv4,Uf_iv4] = make_Hankel(u_iv4,p,f);
[~,Yp_iv4,~] = make_Hankel(y_iv4,p,f);
Ziv4 = [Up_iv4;Yp_iv4;Uf_iv4];

%% 5) approximation of optimal IV w/ controller information
close all;
% checking deadbeat observer for controller

% controller-dependent matrices
[Ac,Bc,Cc,Dc] = ssdata(Cz);
Tf_c = make_blk_tril_toeplitz(Ac,Bc,Cc,Dc,f);
Gf_c = make_ext_obsv(Ac,Cc,f);
Gp_c = make_ext_obsv(Ac,Cc,p);
Kf_c = make_ext_ctrb(Ac,Bc,f,rev=true);
Acf  = Ac^f;

% get predictor form with eigs(Ac-Kc*Cc) \approx 0
nxc = size(Ac,1);
Kc = acker(Ac.',Cc.',zeros(nxc,1)).';
tAc = Ac-Kc*Cc;
tBc = Bc-Kc*Dc;
tKp_ce = make_ext_ctrb(tAc,tBc,p,rev=true);
tKp_cu = make_ext_ctrb(tAc,Kc,p,rev=true);

xc1_est = pinv(Gp_c)*Up_r1...
         -pinv(Gp_c)*make_blk_tril_toeplitz(Ac,Bc,Cc,Dc,p)*(Rp_r0-Yp_r1);
Uf2_r1 = Gf_c*tKp_ce*(Rp_r0-Yp_r1) + Gf_c*tKp_cu*Up_r1 + Tf_c*(Rf_r0-Yf_r1);

% -----------------------------------------

[u_iv5,y_iv5,xCz_iv] = approx_IV_controller_info(u1,y1,y_ref(:,1:Nbar),p,f,Ziv3,Cz,1);
[~,Up_iv5,Uf_iv5] = make_Hankel(u_iv5,p,f);
[~,Yp_iv5,~] = make_Hankel(y_iv5,p,f);
Ziv5 = [Up_iv5;Yp_iv5;Uf_iv5];

% ------------------------------- plotting --------------------------------
figure();
tl = tiledlayout(2,1,'TileSpacing','compact');
ax1_11 = nexttile;
plot(y1); hold on; plot(y2,'LineWidth',2); plot(y_iv4); plot(y_iv5); plot(y_ref); plot(e);
legend({'$y$','$y_{iv}^*$','$y_{iv}^{nc}$','$y_{iv}^{c}$','ref','$e$'},'interpreter','latex');
% ylim([-20 20]);
ax1_12 = nexttile;
plot(u1); hold on; plot(u2,'LineWidth',2); plot(u_iv4); plot(u_iv5);
% ylim([-20 20]);
linkaxes([ax1_11,ax1_12],'x');
% xlim([0 40]);

% --------------------------- get SPC controllers -------------------------
return
% 1) optimal IV
Lf1 = Yf_r1*Zopt.'*pinv(W1*Zopt.');
[Cz1,x0_SPC] = Lf_2_SPC(Lf1,up2usol,yp2usol,urf2usol,yrf2usol,nu,ny,p,f,up=up1,yp=yp1,urf=urf,yrf=yrf);

% 2) SPC using open-loop IV
Lf2 = Yf_r1*W1.'*pinv(W1*W1.');
Cz2 = Lf_2_SPC(Lf2,up2usol,yp2usol,urf2usol,yrf2usol,nu,ny,p,f);

% 3) SPC using LCF
Lf3 = Yf_r1*Ziv3.'*pinv(W1*Ziv3.');
Cz3 = Lf_2_SPC(Lf3,up2usol,yp2usol,urf2usol,yrf2usol,nu,ny,p,f);

% 4) SPC using approximation of optimal IV w/o controller information
Lf4 = Yf_r1*Ziv4.'*pinv(W1*Ziv4.');
Cz4 = Lf_2_SPC(Lf4,up2usol,yp2usol,urf2usol,yrf2usol,nu,ny,p,f);
