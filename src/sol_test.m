%% testing solution to EIV problem with H(q\inv)=F(q\inv) = I_{fl}
close all; clc; clear; rng default;
opts.save = false;
opts.plot = false;
p = 3;
f = 4;
N = 1e5;
Ncl = 1e4;
dRk = 0;
Rk = 1;
Qk = 1e3;
seed = 2;
rng(seed);
% ------------------------- add relevant paths ----------------------------
% go to, save, and add to path src directory
[src_dir, ~  , ~] = fileparts(which(mfilename)); cd(src_dir);
addpath(genpath(src_dir));

% go to, save, and add to path bin directory (and paths to subdirectories)
cd('..'); cd('bin'); bin_dir = pwd; addpath(bin_dir);
addpath(fullfile(bin_dir,'external','casadi-v3.6.7'));   % <- add path for 1) here
if opts.plot
addpath(fullfile(bin_dir,'external','crameri_colours')); % <- add path for 2) here
end

if opts.save
    % go to save, and add to path raw data directory
    cd(src_dir); cd('..'); cd(append('data',filesep,'raw')); data_dir = pwd;
    if isfield(opts,'raw_dir')
        % add opts.raw_dir to raw data path if it exists
        data_dir = fullfile(data_dir,opts.raw_dir{:});
        if ~isfolder(data_dir) % make folder if it doesn't exist yet
            mkdir(data_dir);
        end
    end
    addpath(data_dir);
end

% go back to src directory
cd(src_dir);

%% Simulation settings
fprintf('Setting simulation settings...\n');

opts.sys = 3;
% plant model
switch opts.sys
    case 0
        nu = 3; ny = 2; nx = 4;
        sys = drss(nx,ny,nu);
        [A,B,C,D] = ssdata(sys); 
        K0 = place(A.',C.',linspace(0.8,1,nx)).';
        A = A-K0*C; B = B-K0*D; clear K0;
        K = place(A.',C.',linspace(0.9,0.95,nx)).';
        plant = ss(A,[B,K],C,[D, eye(ny)],1);
        % Cz0 = eye(nu)*pid(1,0.1,0,1e3,plant.Ts);
        % Cz0 = ss(Cz0);
        fn_Cz0 = 'Cz0_random.mat';
        W1 = makeweight(33,5,0.5);  W1 = c2d(W1,plant.Ts,'tustin');
        W3 = makeweight(0.5,20,20); W3 = c2d(W3,plant.Ts,'tustin');
        W2 = [];
    case 1
        [plant,nx,nu,ny,A,B,C,D,K,~] = model_Landau1995();
        W1 = makeweight(33,5,0.5);  W1 = c2d(W1,plant.Ts,'tustin');
        W3 = makeweight(0.5,20,20); W3 = c2d(W3,plant.Ts,'tustin');
        W2 = [];
        fn_Cz0 = 'Cz0_Landau1995.mat';
    case 2
        [plant,nx,nu,ny,A,B,C,D,K,~] = model_Bemporad2002(At_poles=[0.95, 0.9]);
        plant.Ts = 1;
        W1 = makeweight(db2mag(80),[pi/plant.Ts*0.9 1],0.5,plant.Ts);
        W3 = makeweight(0.5,[pi/plant.Ts*0.95 1],20,plant.Ts);
        W2 = ss(1e-1);
        fn_Cz0 = 'Cz0_Bemporad2002.mat';
    case 3
        % system from Favoreel 1999; SPC: Subspace Predictive Control
        A = [ 4.40 1 0 0 0;
             -8.09 0 1 0 0;
              7.83 0 0 1 0;
             -4.00 0 0 0 1;
              0.86 0 0 0 0];
        B = [0.00098 0.01299 0.01859 0.0033 -0.00002].';
        C = eye(1,5);
        D = 0;
        K = [2.3 -6.64 7.515 -4.0146 0.86336].';
        
        nx = size(A,1); nu = size(B,2); ny = size(C,1);
        
        plant = ss(A,[B K], C, [D eye(ny,nu)],1);
        W1 = makeweight(db2mag(80),[pi/plant.Ts*0.58 1],0.85,plant.Ts);
        W3 = makeweight(0.85,[pi/plant.Ts*0.60 1],20,plant.Ts);
        W2 = ss(1e-1);
        fn_Cz0 = 'Cz0_Favoreel1999.mat';
    case 4
        [plant,Cz0,nx,nu,ny,A,B,C,D,K,Re] = model_Wang2023();
        fn_Cz0 = 'Cz0_Wang2023.mat';
        plant.Ts = 1;
        W1 = makeweight(db2mag(80),[pi/plant.Ts*0.58 1],0.85,plant.Ts);
        W3 = makeweight(0.85,[pi/plant.Ts*0.60 1],20,plant.Ts);
        W2 = [];%ss(1e-1);
end

% naming signals
uk_name = arrayfun(@(j) sprintf('u0_%d', j), 1:nu, 'UniformOutput', false);
ek_name = arrayfun(@(j) sprintf('e0_%d', j), 1:ny, 'UniformOutput', false);
yk_name = arrayfun(@(j) sprintf('y0_%d', j), 1:ny, 'UniformOutput', false);
plant.u(1:nu)     =  uk_name;
plant.u(nu+1:end) =  ek_name;
plant.y = yk_name;

% saving to opts structure
[opts.ny,opts.nu,opts.plant] = deal(ny,nu,plant);

Re = tril(rand(ny),-1);
Re = Re+diag(abs(rand(ny,1)+0.2));
Re = sqrtm(Re*Re.'); Re = Re*1e-1;

% ----------------------- make/load initial controller --------------------
switch opts.sys
    case {0,1,2,3,4}
        if isfile(fn_Cz0)
            Cz0 = load(fn_Cz0).Cz0;
        else
            [Cz0,~,~,~] = mixsyn(plant(:,1:nu),W1,[],W3);
        end
end
% naming signals
Cz0.u = arrayfun(@(j) sprintf('er0_%d', j), 1:ny, 'UniformOutput', false);
Cz0.y = plant.u(1:nu);

% ---------------------- set SPC controller settings ----------------------
% weights
dRk= dRk*eye(nu);  dR = kron(speye(f),dRk);
Rk = Rk*eye(nu);   R = kron(speye(f),Rk);
Qk = Qk*eye(ny);   Q = kron(speye(f),Qk);

% ----------------------- set simulation lengths --------------------------
Nbar = p + f + N -1; % sim. length of initial controller

% ------------------------------ set references ---------------------------
yr0  = make_reference(Nbar,ny);
% yr0 = idinput(Nbar,'prbs',[],[-1 1]).'; yr0 = repmat(yr0,ny,1);
[~,~,Rf_yr0] = make_Hankel(yr0,p,f);

% -> references for SPC controllers
yr1 = make_reference(Ncl+f,ny);%,var_dyr=0.01); % y-ref
P0  = dcgain(plant(:,1:nu));                    % DC gain
ur1 = P0\yr1;                                   % u-ref

%% Data generation (closed-loop)
fprintf('Obtaining initial closed-loop data...\n');

% make closed-loop system
fbsum = cell(ny,1);
for k = 1:ny
    fbsum{k} = sumblk(sprintf('er0_%d = r0_%d - y0_%d', k,k,k));
end
rk_name = arrayfun(@(j) sprintf('r0_%d', j), 1:ny, 'UniformOutput', false); % r_k
conOpts = connectOptions("Simplify",false);
Tcl0 = connect(Cz0,plant,fbsum{:},[rk_name(:).',ek_name(:).'],[uk_name(:).',yk_name(:).'],conOpts);

% ---------------------- simulation without noise -------------------------
[uy0_2,~,~] = lsim(Tcl0,[yr0;zeros(ny,Nbar)],[]); uy0_2 = uy0_2.';
u0_2 = uy0_2(1:nu,:); y0_2 = uy0_2(nu+1:end,:); clear uy0_2;

% create Hankel matrices
[~,Up_r02,Uf_r02] = make_Hankel(u0_2,p,f); % - data w/o noise
[~,Yp_r02,Yf_r02] = make_Hankel(y0_2,p,f);

% ----------------------- simulation with noise ---------------------------
e0 = mvnrnd(zeros(ny,1),Re,Nbar).'; % create innovation noise
[uy0,~,xcl0] = lsim(Tcl0,[yr0;e0],[]); uy0 = uy0.'; xcl0 = xcl0.';
u0 = uy0(1:nu,:); y0 = uy0(nu+1:end,:); clear uy0; % get inputs and outputs
xcl0_plant = xcl0(size(Cz0.A,1)+1:size(Cz0.A,1)+nx,end);  % get final state of plant

% create data matrices
[~,Up_r01,Uf_r01] = make_Hankel(u0,p,f); % - data w/ noise
[~,Yp_r01,Yf_r01] = make_Hankel(y0,p,f);
[~,Ep_r01,Ef_r01] = make_Hankel(y0,p,f);


%%
[Up2Yf_inno, Yp2Yf_inno, Uf2Yf_inno] = get_actual_matrices(A,B,C,D,K,p,f);
Lf = [Up2Yf_inno, Yp2Yf_inno, Uf2Yf_inno];
Hf = make_blk_tril_toeplitz(A,K,C,eye(ny),f);

Zopt = [Up_r01;Yp_r01;Uf_r02];
Psi = [Up_r01;Yp_r01;Uf_r01];
Lf_est = Yf_r01*Zopt.'/(Psi*Zopt.');

HfEf_est = Yf_r01-Lf_est*Psi;

% estimating Hf
[~,Hf_est1] = qr(HfEf_est.','econ'); Hf_est1 = Hf_est1.';
Hf_est2 = chol(HfEf_est*HfEf_est.',"lower")
Cov2 = cov((Hf_est2\HfEf_est).')
Hf_est3 = Hf_est2/Cov2
Cov3 = cov((Hf_est3\HfEf_est).')

%%
% make Page matrices
Yf = Yf_r01(:,1:f:end);
Uf = Uf_r01(:,1:f:end);
Yp = Yp_r01(:,1:f:end);
Up = Up_r01(:,1:f:end);
tUf= Uf_r02(:,1:f:end);
Ef = Ef_r01(:,1:f:end);

N = size(Ef,2);




Psi = [Up; Yp; Uf];
tPsi= [Up; Yp; tUf];
model_err = Yf-Lf*Psi-Hf*Ef;

Sigma_eps = Hf*kron(eye(f),Re)*Hf.';
Lf_est1 = get_Lf_est(Sigma_eps,f,ny,tPsi,Psi,Yf);
M = rand(f*ny,f*ny); M = (M+M.')/2;
if min(eig(M))<= 0
    M = M + eye(size(M))*abs(min(eig(M)))*10;
end
S2 = M*Sigma_eps*M;
Lf_est2 = get_Lf_est(S2,f,ny,tPsi,Psi,Yf);
Lf_est3 = Yf*tPsi.'*pinv(Psi*tPsi.');

%
plot(y0.'); hold on; plot(yr0.')

%%
clear Yf Uf Yp Up tUf Ef

N = 1e5;
tUpf = randn(nu*(p+f),N);
Upf = tUpf + randn(nu*(p+f),N)*0.01;
Up  = Upf(1:p*nu,:);
Uf  = Upf(p*nu+1:end,:);
tUf = tUpf(p*nu+1:end,:);
Yp = randn(ny*p,N);
Ef = mvnrnd(zeros(f*ny,1),kron(eye(f),Re),N).';

Psi  = [Up;Yp;Uf];
tPsi = [Up;Yp;tUf];

Hf = eye(4);
Hf(3:4,1:2) = rand(2,2);
Sigma_eps = Hf*kron(eye(f),Re)*Hf.';
inv_Sigma_eps = inv(Sigma_eps);
Yf = Lf*Psi + Hf*Ef;

Lf_est1 = get_Lf_est(Sigma_eps,f,ny,tPsi,Psi,Yf);
Lf_est2 = get_Lf_est(S2,f,ny,tPsi,Psi,Yf);
Lf_est3 = Yf*tPsi.'*pinv(Psi*tPsi.');

%%
a = 0;
b= 0;
kmax = N;
for k = 1:kmax
    a = a + kron(speye(f*ny),tPsi(:,k))*inv_Sigma_eps*Yf(:,k);
    b = b + tPsi(:,k)*Yf(:,k).';
end
a = a/kmax;
b = kron(speye(f*ny),b)*reshape(inv_Sigma_eps,[],1)/kmax;
(a-b)./a*100

%%
function Lf_est = get_Lf_est(S,f,ny,tPsi,Psi,Yf)
N = size(Yf,2);
inv_S = inv(S);
Sum1 = 0;
Sum2 = 0;
for ki = 1:N
    Sum1 = Sum1 + kron(speye(f*ny),tPsi(:,ki))/S*kron(speye(f*ny),Psi(:,ki).');
    Sum2 = Sum2 + kron(speye(f*ny),tPsi(:,ki))/S*Yf(:,ki);
end
Sum1 = Sum1/N;
Sum2 = Sum2/N;
theta_est = Sum1\Sum2;
Lf_est = reshape(theta_est,size(Psi,1),f*ny).';
end