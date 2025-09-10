function main(opts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           DDPC using an Optimal-IV
%           Authors: R. Dinkla, T. Oomen, J.W. van Wingerden
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
arguments (Input)
    opts.Re   (1,1) double  = 4.81e-2;  % innovation noise variance
    opts.plot       logical = false;
    opts.p    (1,1) double  = 20;       % window lengths
    opts.f    (1,1) double  = 20;
    opts.N    (1,1) double  = 1e4;      % number of data matrix columns
    opts.Ncl  (1,1) double  = 1500;     % simulation length of SPC
    opts.dRk  (1,1) double  = 1;        % weights
    opts.Rk   (1,1) double  = 1;
    opts.Qk   (1,1) double  = 1e2;
    opts.seed (1,1) double  = 1;
    opts.save       logical = true;     % save data
    opts.raw_dir    cell;               % subdirectory of raw data directory in which to save files
    opts.sys  (1,1) double = 1;         % flag for model selection
end
[Re, p, f, N, Ncl, dRk, Rk, Qk, seed] = deal(opts.Re, opts.p, opts.f, opts.N, opts.Ncl, opts.dRk, opts.Rk, opts.Qk, opts.seed);
rng(seed);

% Requirements:
% 1) Casadi                                     v3.6.7
% 2) crameri_colours                            v1.09
%    Used for plotting if opts.plot = true. Obtained from
%    https://nl.mathworks.com/matlabcentral/fileexchange/68546-crameri-perceptually-uniform-scientific-colormaps
% 3) Control System Toolbox                     v24.2
% 4) System Identification Toolbox              v24.2
% 5) Robust Control Toolbox                     v24.2
% 6) Statistics and Machine Learning Toolbox    v24.2

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
    % go to, save, and add to path raw data directory
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

% plant model
switch opts.sys
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
[opts.ny,opts.nu,opts.uk_name,opts.yk_name] = deal(ny,nu,uk_name,yk_name);

% ----------------------- make/load initial controller --------------------
if isfile(fn_Cz0)
    Cz0 = load(fn_Cz0).Cz0;
else
    [Cz0,~,~,~] = mixsyn(plant(:,1:nu),W1,[],W3);
    save(fn_Cz0,"Cz0");
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
% yr0 = idinput(Nbar,'prbs',[],[-1 1]).';

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

% ----------------------- simulation with noise ---------------------------
e0 = mvnrnd(zeros(ny,1),Re,Nbar).'; % create innovation noise
[uy0,~,xcl0] = lsim(Tcl0,[yr0;e0],[]); uy0 = uy0.'; xcl0 = xcl0.';
u0 = uy0(1:nu,:); y0 = uy0(nu+1:end,:); clear uy0; % get inputs and outputs
xcl0_plant = xcl0(size(Cz0.A,1)+1:size(Cz0.A,1)+nx,end);  % get final state of plant

% create data matrices
[~,Up_r01,Uf_r01] = make_Hankel(u0,p,f); % - data w/ noise
[~,Yp_r01,Yf_r01] = make_Hankel(y0,p,f);
[~,Ep_r01,Ef_r01] = make_Hankel(y0,p,f);

% ------------------- simulations without future noise --------------------
Uf_r02 = nan(nu*f,N);
Yf_r02 = nan(ny*f,N);
tic
for kN = 1:N
    kk = p+kN;
    uy_f_iv_opt = lsim(Tcl0,[yr0(:,kk:kk+f-1);zeros(ny,f)],[],xcl0(:,kk).').';
    Uf_r02(:,kN) = reshape(uy_f_iv_opt(1:nu,:),nu*f,1);
    Yf_r02(:,kN) = reshape(uy_f_iv_opt(nu+(1:ny),:),ny*f,1);
end
toc

%% Instrumental Variable Matrices
% determine past data length to use - varrho
[Ac,Bc,Cc,Dc] = ssdata(Cz0);
nxc = size(Ac,1);
Gamma_c = make_ext_obsv(Ac,Cc,nxc);
if rank(Gamma_c)~= nxc
    error('Controller is not observable');
end
% rho = lag of the controller
rho = nan;
for k = 1:nxc
    if rank(Gamma_c(1:nu*k,:)) == nxc
        rho = k;
        break;
    end
end
opts.rho = rho; % save to structure

% ============= create data structure and explain cases ===================
%    Name   |       Description
% ----------|----------------------------------------------------------
%   iv1     | SPC using open-loop IV
%   iv2a    | SPC using optimal IV
%   iv2b    | SPC using optimal IV + Yf_iv
%   iv2c    | SPC using optimal IV + Yf_iv + 2SLS
%   iv3a    | SPC using LCF-IV
%   iv3b    | SPC using LCF-IV + 2SLS
%   iv4a    | SPC using approx. opt. IV w/o controller info.
%   iv4b    | SPC using approx. opt. IV w/o controller info. + Yf_iv
%   iv4c    | SPC using approx. opt. IV w/o controller info. + Yf_iv + 2SLS
%   iv5a    | SPC using approx. opt. IV w/  controller info.
%   iv5b    | SPC using approx. opt. IV w/  controller info. + Yf_iv
%   iv5c    | SPC using approx. opt. IV w/  controller info. + Yf_iv + 2SLS
%   iv6a    | SPC using basic IV: future reference
%   iv6b    | SPC using basic IV: future reference + 2SLS
%   CLSPC   | CL-SPC
%   actLf   | SPC using the actual matrix Lf

% make structure array for data
Cases = {'iv1','iv2a','iv2b','iv2c','iv3a','iv3b','iv4a','iv4b','iv4c',...
         'iv5a','iv5b','iv5c','iv6a','iv6b','CLSPC','actLf'};
Descr = {...
'open-loop IV',...                                        iv1    + SPC
'optimal IV',...                                          iv2a   + SPC
'optimal IV + Yf_iv',...                                  iv2b   + SPC
'optimal IV + Yf_iv + 2SLS',...                           iv2c   + SPC
'LCF-IV',...                                              iv3a   + SPC
'LCF-IV + 2SLS',...                                       iv3b   + SPC
'approx. opt. IV w/o controller info.',...                iv4a   + SPC
'approx. opt. IV w/o controller info. + Yf_iv',...        iv4b   + SPC
'approx. opt. IV w/o controller info. + Yf_iv + 2SLS',... iv4c   + SPC
'approx. opt. IV w/  controller info.',...                iv5a   + SPC
'approx. opt. IV w/  controller info. + Yf_iv',...        iv5b   + SPC
'approx. opt. IV w/  controller info. + Yf_iv + 2SLS',... iv5c   + SPC
'basic IV: future reference',...                          iv6a   + SPC   
'basic IV: future reference + 2SLS',...                   iv6a   + SPC
'CL-SPC',...                                              CLSPC
'SPC using the actual matrix Lf'};%                       actLf  + SPC
nCz = numel(Cases);

% ========================== calculate IVs ================================
% some necessary calculations before assigning IVs using IV_4_DDPC class
fprintf('Obtaining IVs...\n');

% ---------------------- (3) LCF-IV ---------------------------------------
% see "Data-Driven Predictive Control Using  Closed-Loop Data: An
% Instrumental  Variable Approach" (2023) by Wang et al.
% DOI: 10.1109/LCSYS.2023.3340444
[~,Vc,Uc] = lncf(ss(Cz0));
Hcv = make_blk_tril_toeplitz(Vc.A,Vc.B,Vc.C,Vc.D,f);
Hcu = make_blk_tril_toeplitz(Uc.A,Uc.B,Uc.C,Uc.D,f);
IV_Theta = Hcv*Uf_r01 + Hcu*Yf_r01;

% ---------------------- (4) w/o controller info. -------------------------
% w/o -> don't known Cz0 exactly, but know rho & feedback configuration
[Uf_iv4,Yf_iv4] = approx_IV_no_controller_info(u0,y0,yr0,opts);

% ---------------------- (5) w/ controller info. --------------------------
[Uf_iv5,Yf_iv5] = approx_IV_controller_info(u0,y0,yr0,opts,Cz0);

% =========================== assign IVs ==================================
% using an instance of the IV_4_DDPC class

% ---------------------- (1) open-loop IV ---------------------------------
% 1) i.e. 'normal' least-squares regression)
Z = IV_4_DDPC(u0,y0,p,f); % initializes IV object & makes open-loop IV ('iv1')

for kIV = 2:nCz-2 % loop over remaining IV names
    IV_name = Cases{kIV};    % IV name
    IV_descr = Descr{kIV};   % IV description

    switch IV_name
% ---------------------- (2) optimal IV -----------------------------------
        case 'iv2a'         % 2a) w/o Yf_iv (result for minimum asymptotic variance)
            Z0 = Uf_r02;
            TSLS = false;
        case 'iv2b'         % 2b) = iv2a + Yf_iv
            Z0 = [Uf_r02;Yf_r02];
            TSLS = false;
        case 'iv2c'         % 2c) = iv2a + Yf_iv + 2SLS
            Z0 = [Uf_r02;Yf_r02];
            TSLS = true;
    
% ---------------------- (3) LCF-IV ---------------------------------------
        case 'iv3a'         % 3a) orignal form of LCF-IV
            Z0 = [IV_Theta; Rf_yr0];
            TSLS = false;
        case 'iv3b'         % 3b) => 3a) + 2SLS
            Z0 = [IV_Theta; Rf_yr0];
            TSLS = true;
    
% ---------------------- (4) w/o controller info. -------------------------
        case 'iv4a'         % 4a) without Yf_iv
            Z0 = Uf_iv4;
            TSLS = false;
        case 'iv4b'         % 4b) => 4a) + Yf_iv
            Z0 = [Uf_iv4;Yf_iv4];
            TSLS = false;
        case 'iv4c'         % 4c) => 4a) + Yf_iv + 2SLS
            Z0 = [Uf_iv4;Yf_iv4];
            TSLS = true;
    
% ---------------------- (5) w/ controller info. --------------------------
        case 'iv5a'         % 5a) without Yf_iv
            Z0 = Uf_iv5;
            TSLS = false;
        case 'iv5b'         % 5b) => 5a) + Yf_iv
            Z0 = [Uf_iv5;Yf_iv5];
            TSLS = false;
        case 'iv5c'         % 5c) => 5a) + Yf_iv + 2SLS
            Z0 = [Uf_iv5;Yf_iv5];
            TSLS = true;
    
% ---------------------- (6) basic IV -------------------------------------
        case 'iv6a'         % 6a) IV composed of future reference signal
            Z0 = Rf_yr0;
            TSLS = false;
        case 'iv6b'         % 6b) => 6a) + SLS
            Z0 = Rf_yr0;
            TSLS = true;
    
    end
    
    % create IV based on Z0, TSLS flag, and description
    Z.add_IV(IV_name,Z0,TSLS=TSLS,descr=IV_descr);

end

% ============= calc. mean & std. dev. of IV trajectories =================
% ---------------------- (1) open-loop IV ---------------------------------
% see u0 & y0, std. dev. = 0

% ---------------------- (2) optimal IV -----------------------------------
[u_iv2a,u_iv2a_std] = diag_stats(Uf_r02,nr=nu,anti=true); % w/o Yf_iv
[y_iv2b,y_iv2b_std] = diag_stats(Yf_r02,nr=ny,anti=true); % w/  Yf_iv
Uf_iv2c = Z.iv2c_;                                        % w/  Yf_iv + 2SLS
[u_iv2c,u_iv2c_std] = diag_stats(Uf_iv2c,nr=nu,anti=true); 

% ---------------------- (3) LCF-IV ---------------------------------------
[y_iv3a,y_iv3a_std] = diag_stats(IV_Theta,nr=nu,anti=true); % interpretation of IV_Theta
Uf_iv3b = Z.iv3b_;                                        % w/ 2SLS
[u_iv3b,u_iv3b_std] = diag_stats(Uf_iv3b,nr=nu,anti=true);

% ---------------------- (4) w/o controller info. -------------------------
[u_iv4a,u_iv4a_std] = diag_stats(Uf_iv4,nr=nu,anti=true); % w/o Yf_iv
[y_iv4b,y_iv4b_std] = diag_stats(Yf_iv4,nr=ny,anti=true); % w/  Yf_iv
Uf_iv4c = Z.iv4c_;                                        % w/  Yf_iv + 2SLS
[u_iv4c,u_iv4c_std] = diag_stats(Uf_iv4c,nr=nu,anti=true);

% ---------------------- (5) w/ controller info. --------------------------
[u_iv5a,u_iv5a_std] = diag_stats(Uf_iv5,nr=nu,anti=true); % w/o Yf_iv
[y_iv5b,y_iv5b_std] = diag_stats(Yf_iv5,nr=ny,anti=true); % w/  Yf_iv
Uf_iv5c = Z.iv5c_;                                        % w/  Yf_iv + 2SLS
[u_iv5c,u_iv5c_std] = diag_stats(Uf_iv5c,nr=nu,anti=true);

% ---------------------- (6) basic IV --------------------------
y_iv6a = yr0(:,p+1:end);                                  % reference
Uf_iv6b = Z.iv6b_;                                        % reference + 2SLS
[u_iv6b,u_iv6b_std] = diag_stats(Uf_iv6b,nr=nu,anti=true);

%% Get Subspace Predictive Controllers
% ---------------------- get (CL-)SPC controllers -------------------------
fprintf('Obtaining controllers...\n');

% make functions to get SPC controllers
usol_funs = get_solver(nu,ny,p,f,Q,R,dR); % structure w/ up2usol, yp2usol, urf2usol, yrf2usol

% name channels in SPC controllers
urkf_name = arrayfun(@(j) sprintf('ur%d_%d',f, j), 1:nu, 'UniformOutput', false);
yrkf_name = arrayfun(@(j) sprintf('yr%d_%d',f, j), 1:ny, 'UniformOutput', false);
[opts.urkf_name,opts.yrkf_name] = deal(urkf_name,yrkf_name);

% for initial state of SPC controllers
up0 = u0(:,end-p+1:end); up0 = up0(:);  % past input data
yp0 = y0(:,end-p+1:end); yp0 = yp0(:);  % past output data
urf = ur1(:,1:f);        urf = urf(:);  % future input references
yrf = yr1(:,1:f);        yrf = yrf(:);  % future output references

% create controllers
for iCz = 1:nCz
    Czn = Cases{iCz}; % name of controller
    switch Czn
% -------------------------- (1-6) SPCs based on an IV --------------------
        case Cases(1:nCz-2)
            Lf.(Czn) = Yf_r01*Z.(Czn).'*pinv(Z.iv1*Z.(Czn).');
            if iCz > 1
                Cz.(Czn) = Lf_2_SPC(Lf.(Czn),usol_funs,opts);
            else
                % also create initial state of the controller
                [Cz.(Czn),x1_0_SPC] = Lf_2_SPC(Lf.(Czn),usol_funs,opts,up=up0,yp=yp0,urf=urf,yrf=yrf);
            end

% -------------------------- (7) CL-SPC -----------------------------------
        case 'CLSPC'
            Lf.(Czn) = get_Lf_CL_SPC(u0,y0,p,f,nu,ny);
            Cz.(Czn) = Lf_2_SPC(Lf.(Czn),usol_funs,opts);

% -------------------------- (8) SPC w/ actual Lf -------------------------
        case 'actLf'
            [Up2Yf_inno, Yp2Yf_inno, Uf2Yf_inno] = get_actual_matrices(A,B,C,D,K,p,f);
            Lf.(Czn) = [Up2Yf_inno Yp2Yf_inno Uf2Yf_inno];
            Cz.(Czn) = Lf_2_SPC(Lf.(Czn),usol_funs,opts);
    end
end

%% Run closed-loop simulations
fprintf('Running closed-loop simulations...\n');

e1 = mvnrnd(zeros(ny,1),Re,Ncl).'; % create innovation noise
x1_0_plant = plant.A*xcl0_plant+plant.B*[u0(:,end);e0(:,end)];
x1_0_cl = [x1_0_SPC;x1_0_plant];

% create closed-loop systems
Tcl_in  = [urkf_name(:)',yrkf_name(:)',ek_name(:)'];
Tcl_out = [uk_name(:)',yk_name(:)'];
Tcls = cell(nCz,1);
for kCz = 1:nCz
    Czn = Cases{kCz};
    Tcl.(Czn) = connect(Cz.(Czn),plant,Tcl_in,Tcl_out,conOpts);
end

% define data structures
for kCz = 1:nCz
    Czn = Cases{kCz}; % controller name
    uy_cl = lsim(Tcl.(Czn),[ur1(:,f+1:end);yr1(:,f+1:end);e1],[],x1_0_cl).';
    u_cl.(Czn) = uy_cl(1:nu,:);
    y_cl.(Czn) = uy_cl(nu+1:end,:);
end
fprintf('Closed-loop simulations finished!\n');

%% Calculate average stage costs
fprintf('Calculate average stage costs...\n');

% -------------------------- create functions -----------------------------
% average of (u_k-ur_k).' * Rk * (u_k-ur_k) over Ncl steps
cost_fun_uk = @(u) reshape(ur1(:,1:Ncl)-u,[],1).'*kron(speye(Ncl),Rk)*reshape(ur1(:,1:Ncl)-u,[],1);

% average of du_k.' * dRk * du_k over Ncl steps
cost_fun_duk = @(u) reshape(u-[u0(:,end) u(:,2:end)],[],1).'*kron(speye(Ncl),dRk)*reshape(u-[u0(:,end) u(:,2:end)],[],1);

% average of (y_k-yr_k).' * Rk * (y_k-yr_k) over Ncl steps
cost_fun_yk = @(y) reshape(yr1(:,1:Ncl)-y,[],1).'*kron(speye(Ncl),Qk)*reshape(yr1(:,1:Ncl)-y,[],1);

% -------------------------- perform calculations -------------------------
for kCz = 1:nCz
    Czn = Cases{kCz};
    
    cost_u1.(Czn) = cost_fun_uk( u_cl.(Czn))/Ncl;
    cost_u2.(Czn) = cost_fun_duk(u_cl.(Czn))/Ncl;
    cost_y.(Czn)  = cost_fun_yk( y_cl.(Czn))/Ncl;
    cost_u.(Czn)  = cost_u1.(Czn) + cost_u2.(Czn);
    cost_tot.(Czn)= cost_u.(Czn)  + cost_y.(Czn);  % total cost
end

%% Identification error analysis
fprintf('Calculate identification errors...\n');

FroIDerror = zeros(kCz,3);
for kCz = 1:nCz
    Czn = Cases{kCz};

    IDerror = Lf.(Czn)-Lf.actLf;
    cols = 1:p*nu;
    FroIDerror(kCz,1) = norm(IDerror(:,cols),'fro');
    cols = p*nu+(1:p*ny);
    FroIDerror(kCz,2) = norm(IDerror(:,cols),'fro');
    cols = p*(nu+ny)+(1:f*nu);
    FroIDerror(kCz,3) = norm(IDerror(:,cols),'fro');
end

%% Saving data
if opts.save
    fn = sprintf('seed_%d.mat',seed);
    if ~isfield(opts,'raw_dir')
        % Helper function to trim to minimal digits in scientific notation
        trimmed_exp = @(x) regexprep(sprintf('%e', x), '(\.\d*?)0+(e[+-]?\d+)', '$1$2'); % trims trailing 0s
        % Also remove . if nothing follows
        trimmed_exp = @(x) regexprep(trimmed_exp(x), '\.(e)', '$1');
        
        % Apply formatting
        Re_str  = trimmed_exp(Re); N_str   = trimmed_exp(N);  Ncl_str = trimmed_exp(Ncl);
        Qk_str  = trimmed_exp(Qk); Rk_str  = trimmed_exp(Rk); dRk_str = trimmed_exp(dRk);
        data_dir2 = sprintf('Re_%s_p_%d_f_%d_N_%s_Ncl_%s_Qk_%s_Rk_%s_dRk_%s',Re_str,p,f,N_str,Ncl_str,Qk_str,Rk_str,dRk_str);
        data_dir2 = replace(data_dir2,'.','p');
        data_dir2 = replace(data_dir2,'+','');
        data_dir = fullfile(data_dir,data_dir2);
        if ~isfolder(data_dir)
            % make data folder to store data for different seeds
            mkdir(data_dir);
    
            % also add .txt version of executed main.m file
            copyfile(fullfile(src_dir,append(mfilename,'.m')),...
                     fullfile(data_dir,append(mfilename,'.txt')));
        end
    end
    fn = fullfile(data_dir,fn);
    fprintf('Saving data to file: \n\t%s \n',fn);
    save(fn,'Re','p','f','N','Ncl','Qk','Rk','dRk','seed','Nbar',...
        'plant','Cz0','Tcl0','yr0','yr1','ur1','e0','u0','y0','xcl0',...
        'u0_2','y0_2','Lfs','Czs','Tcls','u_cl','y_cl',...
        'u_iv4a','y_iv4b','u_iv5a','y_iv5b','u_iv9','y_iv9',...
        'cost_u1','cost_u2','cost_u','cost_y','cost_tot',...
        'FroIDerror');
    fprintf('File saved successfully!\n');
end

%% Plotting
if opts.plot
    close all;
% ============================= IV trajectories ===========================
    
    cCram = crameri('roma', 8); % colors

    % free (noisy) response
    [y0_free,~,~] = lsim(plant,[zeros(nu,Nbar);e0],[]); y0_free = y0_free.';
    
    figure(1);
    tiledlayout(2,1,'TileSpacing','compact');
    
    ax1_11 = nexttile; % for Yf_iv
    plot(y_iv6a,'Color','b','LineWidth',2); hold on;                       % reference 
    plot(y0(:,p+1:end),'LineWidth',2,'Color','k');                         % actual output
    plot(e0(:,p+1:end),'Color','r','LineStyle','--');                      % innovation noise
    plotMeanWithStd(y_iv2b,y_iv2b_std,Color=cCram(1,:),LineStyle='-');     % optimal IV
    plotMeanWithStd(y_iv4b,y_iv4b_std,Color=cCram(4,:),LineStyle='--');    % w/o controller info  
    plotMeanWithStd(y_iv5b,y_iv5b_std,Color=cCram(6,:),LineStyle='-');     % w/  controller info
    plotMeanWithStd(y_iv3a,y_iv3a_std,Color=cCram(8,:),LineStyle='-');     % LCF-IV - IV_Theta)
    yLim1 = ax1_11.YLim;
    plot(y0_free(:,p+1:end),'g');                                          % free response
    ylim(yLim1);
    legend('3a,6a) ref',     'actual output',   'noise',...
           '2b) opt. IV', '4b) w/o Cz info', '5b) w/  Cz info',...
           '3a)  LCF-IV Theta','free resp.')
    % legend({'$y$','$y_{iv}^*$','$y_{iv}^{nc}$','$y_{iv}^{c,lcf}$','$y_{iv}^{c,b}$','ref','$e$','$y_{\mathrm{free}}$'},'Interpreter','latex');
    grid on;
    ylabel('$y_k$','Interpreter','latex');
    
    ax1_12 = nexttile; % for Uf_iv
    plotMeanWithStd(u_iv6b,u_iv6b_std,Color='b',lineWidth=2); hold on;     % reference + 2SLS
    plot(u0(:,p+1:end),'LineWidth',2,'Color','k');                         % actual output
    plotMeanWithStd(u_iv2a,u_iv2a_std,Color=cCram(1,:),LineStyle='-');     % optimal IV
    plotMeanWithStd(u_iv2c,u_iv2c_std,Color=cCram(2,:),LineStyle='--');    % optimal IV + Yf_iv + 2SLS
    plotMeanWithStd(u_iv3b,u_iv3b_std,Color=cCram(3,:),LineStyle='-');     % LCF-IV + 2SLS
    plotMeanWithStd(u_iv4a,u_iv4a_std,Color=cCram(4,:),LineStyle='--');    % w/o controller info
    plotMeanWithStd(u_iv4c,u_iv4c_std,Color=cCram(5,:),LineStyle='-');     % w/o controller info + Yf_iv + 2SLS
    plotMeanWithStd(u_iv5a,u_iv5a_std,Color=cCram(6,:),LineStyle='-');     % w/  controller info
    plotMeanWithStd(u_iv5c,u_iv5c_std,Color=cCram(7,:),LineStyle='--');    % w/  controller info + Yf_iv + 2SLS
    legend('6b)  ref + 2SLS',        'actual input',    '2ab) opt. IV',...
           '2c)  opt. IV + 2SLS',    '3b)  LCF + 2SLS', '4ab) w/o Cz info',...
           '4c)  w/o Cz info + 2SLS','5ab) w/ Cz info', '5c)  w/ Cz info + 2SLS')
    grid on;
    ylabel('$u_k$','Interpreter','latex');
    xlabel('Time','Interpreter','latex');
    
    linkaxes([ax1_11,ax1_12],'x');

% ======================= closed-loop simulation data =====================
    % Load a Crameri colormap (e.g., 'batlow')
    colors = crameri('batlow', nCz);  % nCz colors for controllers
    
    % Define line and marker styles
    lineStyles = {'-','--',':','-.'};    % line styles
    markerStyles = {'o','s','^','d','none'};    % include 'none' to skip markers
    
    figure(2);
    tiledlayout(2,1,'TileSpacing','compact');
    
    ax2_1 = nexttile;
    stairs(yr1, 'k-','LineWidth',2,'DisplayName','ref'); hold on;
    
    for kCz = 1:nCz
        Czn = Cases{kCz}; % controller name
        ls = lineStyles{mod(kCz-1,length(lineStyles))+1};
        ms = markerStyles{mod(kCz-1,length(markerStyles))+1};
        stairs(y_cl.(Czn), 'Color', colors(kCz,:), 'LineStyle', ls, ...
             'Marker', ms, 'LineWidth',1.5, 'DisplayName', Czn);
    end
    
    ylim([-15 15]); grid on;
    ylabel('$y_k$','Interpreter','latex');
    legend show
    
    ax2_2 = nexttile;
    stairs(ur1, 'k-','LineWidth',2,'DisplayName','ref'); hold on;
    
    for kCz = 1:nCz
        Czn = Cases{kCz}; % controller name
        ls = lineStyles{mod(kCz-1,length(lineStyles))+1};
        ms = markerStyles{mod(kCz-1,length(markerStyles))+1};
        stairs(u_cl.(Czn), 'Color', colors(kCz,:), 'LineStyle', ls, ...
             'Marker', ms, 'LineWidth',1.5, 'DisplayName', Czn);
    end
    
    ylim([-30 30]); grid on;
    ylabel('$u_k$','Interpreter','latex');
    xlabel('Time', 'Interpreter','latex');
    legend show
    
    linkaxes([ax2_1 ax2_2], 'x');

    
% ======================= visualize identification results ================
    cabsmax = 0;
    for kCz = 1:nCz
        Czn = Cases{kCz}; % controller name
        cabsmax = max(max(abs(Lf.(Czn)),[],"all"),cabsmax);
    end
    
    figure(3);
    tiledlayout(nCz,1,'TileSpacing','compact');
    for kCz = 1:nCz
        Czn = Cases{kCz}; % controller name
        nexttile;
        imagesc_vik(Lf.(Czn),cmax=cabsmax);
    end
    
% ======================= visualize identification error ==================
    figure(4);
    imagesc(FroIDerror);
    cmap = crameri('-davos');
    colormap(cmap);
    colorbar;

% ======================= visualize average stage costs ===================
    figure(5);
    xleg = categorical(Cases);
    xleg = reordercats(xleg,Cases);
    bar(xleg,[cost_u cost_y],'stacked');
    legend({'$\mathcal{J}_u$','$\mathcal{J}_y$'},'interpreter','latex');
    ax_bar = gca;
    ax_bar.TickLabelInterpreter = 'latex';
    max_y = max(cost_tot(~isoutlier(cost_tot)),[],'all');
    ylim(ax_bar,[0 1.05*max_y]);
end
end

%% Helper functions
function Lf = get_Lf_CL_SPC(u1,y1,p,f,nu,ny)
    [~,Up,Uf] = make_Hankel(u1,p,1);
    [~,Yp,Yf] = make_Hankel(y1,p,1);
    L1 = Yf*pinv([Up;Yp;Uf(1:nu,:)]);
    
    C_tKpu_hat = L1(:,1:p*nu);
    C_tKpy_hat = L1(:,p*nu+1:p*(nu+ny));
    D_hat = L1(:,end-nu+1:end);
    
    % construct multiple-step-ahead predictor
    tLest_u = zeros(ny*f,nu*(p+f));
    tLest_u(1:ny,1:nu*(p+1)) = [C_tKpu_hat D_hat];
    tLest_y = zeros(ny*f,ny*(p+f));
    tLest_y(1:ny,1:ny*p) = C_tKpy_hat;
    for kr = 2:f
        tLest_u((kr-1)*ny+1:kr*ny,:) = circshift(tLest_u((kr-2)*ny+1:(kr-1)*ny,:),nu,2);
        tLest_y((kr-1)*ny+1:kr*ny,:) = circshift(tLest_y((kr-2)*ny+1:(kr-1)*ny,:),ny,2);
    end
    tHf = eye(ny*f)-tLest_y(:,end-ny*f+1:end);
    Lf = tHf\[tLest_u(:,1:p*nu) tLest_y(:,1:p*ny) tLest_u(:,end-nu*f+1:end)];
end

function [hFill,hLine] = plotMeanWithStd(mean_vals, std_devs, opts)
% plotMeanWithStd - Plot mean values with shaded standard deviation bands
%
% Syntax:
%   plotMeanWithStd(mean_vals, std_devs)
%   plotMeanWithStd(mean_vals, std_devs, opts)
%
% Description:
%   This function visualizes a mean curve with its uncertainty (standard
%   deviation) as a shaded region. The shaded patch is excluded from the
%   legend, ensuring the legend only describes the mean line. Plot
%   appearance can be customized via the options structure `opts`.
%
% Inputs:
%   mean_vals (1,:) double
%       Vector of mean values to plot.
%
%   std_devs (1,:) double
%       Vector of standard deviations corresponding to `mean_vals`.
%       Must be the same length as `mean_vals`.
%
%   opts (optional) struct with fields:
%       .Color  - Line/patch color. Can be:
%                       * MATLAB short color char (e.g. 'b','r','g')
%                       * RGB triplet (e.g. [0 0.447 0.741])
%                     Default: 'b'
%       .LineStyle - Line style for mean curve (e.g. '-', '--', ':')
%                     Default: '-'
%       .faceAlpha - Transparency of shaded area (0 = transparent,
%                     1 = opaque). Default: 0.2
%       .lineWidth - Width of mean line. Default: 2
%       .x         - x-axis values corresponding to `mean_vals`.
%                     Must be same length as `mean_vals`.
%                     Default: 1:length(mean_vals)
%
% Outputs:
%   hFill - handle to shaded patch (excluded from legend)
%   hLine - handle to mean line (included in legend)
%
% Example:
%   x = linspace(0,2*pi,100);
%   y = sin(x);
%   s = 0.1*ones(size(y));
%   opts.Color = [0.2 0.6 0.8];
%   opts.LineStyle = '--';
%   opts.lineWidth = 3;
%   plotMeanWithStd(y, s, opts);

    arguments
        mean_vals (1,:) double
        std_devs  (1,:) double {mustBeEqualLength(mean_vals,std_devs)}
        opts.Color  = 'b'
        opts.LineStyle = '-'
        opts.faceAlpha (1,1) double {mustBeGreaterThanOrEqual(opts.faceAlpha,0),mustBeLessThanOrEqual(opts.faceAlpha,1)} = 0.2
        opts.lineWidth (1,1) double {mustBePositive} = 2
        opts.x (1,:) double {mustBeEqualLength(opts.x,std_devs)} = 1:length(mean_vals)
    end

    x = opts.x;

    % Upper and lower bounds
    upper = mean_vals + std_devs;
    lower = mean_vals - std_devs;

    % Convert to RGB if short color char is given
    if ischar(opts.Color) || isstring(opts.Color)
        colorRGB = getColorFromChar(opts.Color);
    else
        colorRGB = opts.Color;
    end

    % Fill area between bounds (not shown in legend)
    hFill = fill([x fliplr(x)], [upper fliplr(lower)], colorRGB, ...
        'FaceAlpha', opts.faceAlpha, 'EdgeColor', 'none', ...
        'HandleVisibility','off'); 
    hold on;

    % Plot mean line (included in legend)
    hLine = plot(x, mean_vals, 'Color', colorRGB, ...
        'LineWidth', opts.lineWidth, 'LineStyle', opts.LineStyle);
end

function rgb = getColorFromChar(c)
% getColorFromChar - Convert MATLAB short color names to RGB triplets
    switch char(c)
        case 'y', rgb = [1 1 0];
        case 'm', rgb = [1 0 1];
        case 'c', rgb = [0 1 1];
        case 'r', rgb = [1 0 0];
        case 'g', rgb = [0 1 0];
        case 'b', rgb = [0 0 1];
        case 'w', rgb = [1 1 1];
        case 'k', rgb = [0 0 0];
        otherwise
            error('Unknown color specifier: %s', c);
    end
end

function mustBeEqualLength(a,b)
% mustBeEqualLength - Validation function for equal-length vectors
    if length(a) ~= length(b)
        error('Inputs must be the same length.');
    end
end
