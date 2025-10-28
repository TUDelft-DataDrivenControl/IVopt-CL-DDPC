%% script to check whether the actual Lf matrix also corresponds with
% actual simulated data
clear;
close all;
clc;

%% get plant
[plant,nu,ny,Cz0,Tcl0,opts,sigs] = init_sims(struct('sys',1));
[A,B,C,D,K] = plant2ABCDK(plant);
nx = size(A,1);

%% DDPC settings
p = 5;
f = 20;

% test below
[Lf_Up, Lf_Yp, Lf_Uf] = get_actual_matrices(plant,p,f);
Lf1 = [Lf_Up, Lf_Yp, Lf_Uf];

% matrix for effect of future noise
Hf = make_blk_tril_toeplitz(A,K,C,1,f);

% matrix for effect of past noise
Gp = make_ext_obsv(A,C,p);
Gf = make_ext_obsv(A,C,f);
At = A-K*C;
Hp = make_blk_tril_toeplitz(A,K,C,1,p);
MatEp = -Gf*At^p*pinv(Gp)*Hp;

%% simulate
N = 1e4;
Nbar = N+f+p-1;
r0 = mvnrnd(zeros(ny,1),eye(ny),Nbar).';
e0 = mvnrnd(zeros(nu,1),eye(nu),Nbar).'; e0=e0*1e-1;
% for initial state
up = zeros(nu,p);
yp = zeros(ny,p);
ep = zeros(ny,p);

% simulate actual system (lsim)
uy0  = lsim(Tcl0,[r0;e0].',[]).'; % with noise
u0   = uy0(1:nu,:);
y0   = uy0(nu+1:end,:); 
% u0 = mvnrnd(zeros(nu,1),eye(nu),Nbar).';
% y0 = lsim(plant,[u0;e0].',[]).';

eAll = [ep e0];
uAll = [up u0];
yAll = [yp y0];
[~,Up0,Uf0] = make_Hankel(uAll,p,f);
[~,Ep0,Ef0] = make_Hankel(eAll,p,f);
[~,Yp0,Yf0] = make_Hankel(yAll,p,f);

% simulating w/ matrix Lf1
Yf1_pred = Lf1*[Up0;Yp0;Uf0];        % output predictions
Yf1 = Yf1_pred + MatEp*Ep0 + Hf*Ef0; % prediction + unknown noise contributions

% check that Yf0 & Yf1 are equal:
fprintf('max(abs(Yf0-Yf1))=%.2e\n',max(abs(Yf0-Yf1),[],'all'));

%% Exploring closed-loop Lf estimation

Z0 = [Up0;Yp0;Uf0];
pinvZ = pinv(Z0);

% estimated matrix (no IV)
Lf1_est1 = Yf1_pred*pinvZ; % w/o Ep & Ef
Lf1_est2 = Yf1*pinvZ;      % w/  Ep & Ef
Ep2Lf1_est2 = MatEp*Ep0*pinvZ; % influence of Ep
Ef2Lf1_est2 = Hf*Ef0*pinvZ;    % influence of Ef
cmax = max(abs([Lf1,Lf1_est1,Lf1_est2,Ep2Lf1_est2,Ef2Lf1_est2]),[],'all');

figure();
subplot(5,1,1); imagesc_vik(Lf1,cmax=cmax);      % actual matrix
subplot(5,1,2); imagesc_vik(Lf1_est1,cmax=cmax); % est. w/o Ep & Ef
subplot(5,1,3); imagesc_vik(Lf1_est2,cmax=cmax); % est. w/  Ep & Ef
subplot(5,1,4); imagesc_vik(Ep2Lf1_est2,cmax=cmax); % influence of Ep
subplot(5,1,5); imagesc_vik(Ef2Lf1_est2,cmax=cmax); % influence of Ef

% check predictions
Yf1_pred1 = Lf1_est1*Z0;
Yf1_pred2 = Lf1_est2*Z0;
figure();
subplot(2,1,1); imagesc_vik(Yf1_pred1-Yf1_pred);
subplot(2,1,2); imagesc_vik(Yf1_pred2-Yf1_pred);


% figure()
% subplot(2,1,1); imagesc_vik(Yp0*Ep0.'/N);
% subplot(2,1,2); imagesc_vik(Up0*Ep0.'/N);

