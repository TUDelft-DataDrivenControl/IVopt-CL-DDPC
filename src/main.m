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
    opts.raw_dir    cell;
    opts.sys  (1,1) double = 1;         % flag for model selection
    opts.SNR  (1,1) double              % noise to signal ratio [-] (not in dB)
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
[opts.ny,opts.nu,opts.plant] = deal(ny,nu,plant);

% ----------------------- make/load initial controller --------------------
switch opts.sys
    case {1,2,3,4}
        if isfile(fn_Cz0)
            Cz0 = load(fn_Cz0).Cz0;
        else
            [Cz0,~,~,~] = mixsyn(plant(:,1:nu),W1,[],W3);
            save(fn_Cz0,"Cz0");
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
% yr0 = idinput(Nbar,'prbs',[],[-1 1]).';

wr = lsim(Cz0,yr0).'; % for the form u_k = w_k - Cz(q) y_k
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

% adjust Re to match SNR (if specified)
if isfield(opts,'SNR') && ~isempty(opts.SNR)
    % SNR calculated for OL: input vs. noise
    yu_e = lsim(plant(:,1:nu)*Tcl0(1:nu,ny+1:end),e0, []);
    ye   = lsim(plant(:,nu+1:end),e0, []);
    Pyur = sum(y0_2.^2);
    Pyue = sum(yu_e.^2);
    Pye  = sum(ye.^2);
    SNR = (Pyur + Pyue)/Pye;
    SNR_min = Pyue/Pye; % obtained when Re = Inf
    SNR_lb = SNR_min*1.01;
    if opts.SNR <= SNR_lb
        opts.SNR = SNR_lb;
    end
    SNR_ratio = opts.SNR/SNR;%Pyur/(Pye*opts.SNR-Pyue); % = (target SNR) / (actual SNR)
    Re = Re*SNR_ratio;
    e0 = e0*sqrt(SNR_ratio);
end
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

%% Get Subspace Predictive Controllers
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

% get max of past controller & plant windows
varrho = max(rho,p);

% ----------------------------- get IVs -----------------------------------
fprintf('Obtaining IVs...\n');

% 1) open-loop IV (i.e. least-squares regression)
Zol  = [Up_r01;Yp_r01;Uf_r01];

% 2) optimal IV
% Zopt = [Up_r01;Yp_r01;Uf_r02];
Zopt  = [Up_r01;Yp_r01;Uf_r02;Yf_r02];
Zopt_2= new_IV(Zopt,Zol,p*(nu+ny));
[u_iv,u_iv_std] = diag_stats(flipud(Uf_r02));
[y_iv,y_iv_std] = diag_stats(flipud(Yf_r02));

% 3) composed IV using LCF
[~,Vc,Uc] = lncf(ss(Cz0));
Hcv = make_blk_tril_toeplitz(Vc.A,Vc.B,Vc.C,Vc.D,f);
Hcu = make_blk_tril_toeplitz(Uc.A,Uc.B,Uc.C,Uc.D,f);
IV_Theta = Hcv*Uf_r01 + Hcu*Yf_r01;
Ziv3 = [Up_r01; Yp_r01; IV_Theta; Rf_yr0];
Ziv3_2 = new_IV(Ziv3,Zol,p*(nu+ny));

% 4) approx. optimal IV w/o controller information
[Uf_iv4,Yf_iv4] = approx_IV_no_controller_info(u0,y0,yr0,rho,opts);
% [~,Up_iv4,Uf_iv4] = make_Hankel(u_iv4,p,f);
% [~,Yp_iv4,Yf_iv4] = make_Hankel(y_iv4,p,f);
Ziv4 = [Up_r01;Yp_r01;Uf_iv4;Yf_iv4];
Ziv4_2 = new_IV(Ziv4,Zol,p*(nu+ny));
[u_iv4,u_iv4_std] = diag_stats(flipud(Uf_iv4));
[y_iv4,y_iv4_std] = diag_stats(flipud(Yf_iv4));

% 5) approx. optimal IV w/ controller information - init. LCF
[Uf_iv5,Yf_iv5,xCz_iv5] = approx_IV_controller_info_v2(u0,y0,yr0,rho,p,f,Ziv3,Cz0,1,plant,e0);
Ziv5 = [Up_r01;Yp_r01;Uf_iv5;Yf_iv5];
Ziv5_2 = new_IV(Ziv5,Zol,p*(nu+ny));
[u_iv5,u_iv5_std] = diag_stats(flipud(Uf_iv5));
[y_iv5,y_iv5_std] = diag_stats(flipud(Yf_iv5));

% 8) basic IV
Ziv8 = [Up_r01; Yp_r01; Rf_yr0];
Ziv8_2 = new_IV(Ziv8,Zol,p*(nu+ny));

% 9) approx optimal IV w/ controller information - init. basic IV
[u_iv9,y_iv9,xCz_iv9] = approx_IV_controller_info(u0,y0,yr0,p,f,Ziv8,Cz0,1);
[~,Up_iv9,Uf_iv9] = make_Hankel(u_iv9,p,f);
[~,Yp_iv9,Yf_iv9] = make_Hankel(y_iv9,p,f);
Ziv9 = [Up_r01;Yp_r01;Uf_iv9;Yf_iv9];
Ziv9_2 = new_IV(Ziv9,Zol,p*(nu+ny));

% ---------------------- get (CL-)SPC controllers -------------------------
fprintf('Obtaining controllers...\n');
[up2usol, yp2usol, urf2usol, yrf2usol] = get_solver(nu,ny,p,f,Q,R,dR);

% for initial state of SPC controllers
up0 = u0(:,end-p+1:end); up0 = up0(:);  % past input data
yp0 = y0(:,end-p+1:end); yp0 = yp0(:);  % past output data
urf = ur1(:,1:f); urf = urf(:);         % future input references
yrf = yr1(:,1:f); yrf = yrf(:);         % future output references

% 1) SPC using open-loop IV
Lf1 = Yf_r01*Zol.'*pinv(Zol*Zol.');
Cz1 = Lf_2_SPC(Lf1,up2usol,yp2usol,urf2usol,yrf2usol,nu,ny,p,f);
urkf_name = arrayfun(@(j) sprintf('ur%d_%d',f, j), 1:nu, 'UniformOutput', false);
yrkf_name = arrayfun(@(j) sprintf('yr%d_%d',f, j), 1:ny, 'UniformOutput', false);
Cz1.u(1:ny)          = yk_name;   % y_k
Cz1.u(ny+1:ny+nu)    = urkf_name; % ur_{k+f}
Cz1.u(end-ny+(1:ny)) = yrkf_name; % yr_{k+f}
Cz1.y = uk_name;                  % u_k

% 2) optimal IV
% Lf2 = IV_GLS(Yf_r01,Zopt,Zol);
Lf_3 = IV_GLS_vec(Yf_r01,Zopt,Zol,ny,0);
% Lf2 = Yf_r01*Zopt.'*pinv(Zol*Zopt.');
Lf2_2 = Yf_r01*Zopt_2.'/(Zopt_2*Zopt_2.');
Lf2 = IV_GMM_vec(Yf_r01,Zopt,Zol,ny,nu,p,f);
[Cz2,x1_0_SPC] = Lf_2_SPC(Lf2,up2usol,yp2usol,urf2usol,yrf2usol,nu,ny,p,f,up=up0,yp=yp0,urf=urf,yrf=yrf);
Cz2.u = Cz1.u; Cz2.y = Cz1.y;

% 3) SPC using LCF
Lf3 = Yf_r01*Ziv3.'*pinv(Zol*Ziv3.');
% Lf3 = Yf_r01*Ziv3_2.'/(Ziv3_2*Ziv3_2.');
Cz3 = Lf_2_SPC(Lf3,up2usol,yp2usol,urf2usol,yrf2usol,nu,ny,p,f);
Cz3.u = Cz1.u; Cz3.y = Cz1.y;

% 4) SPC using approximation of optimal IV w/o controller information
% Lf4 = Yf_r01*Ziv4.'*pinv(Zol*Ziv4.');
Lf4 = Yf_r01*Ziv4_2.'/(Ziv4_2*Ziv4_2.');
Cz4 = Lf_2_SPC(Lf4,up2usol,yp2usol,urf2usol,yrf2usol,nu,ny,p,f);
Cz4.u = Cz1.u; Cz4.y = Cz1.y;

% 5) SPC using approximation of optimal IV w/ controller info. - init LCF
% Lf5 = Yf_r01*Ziv5.'*pinv(Zol*Ziv5.');
Lf5 = Yf_r01*Ziv5_2.'/(Ziv5_2*Ziv5_2.');
Cz5 = Lf_2_SPC(Lf5,up2usol,yp2usol,urf2usol,yrf2usol,nu,ny,p,f);
Cz5.u = Cz1.u; Cz5.y = Cz1.y;

% 6) CL-SPC
Lf6 = get_Lf_CL_SPC(u0,y0,p,f,nu,ny);
Cz6 = Lf_2_SPC(Lf6,up2usol,yp2usol,urf2usol,yrf2usol,nu,ny,p,f);
Cz6.u = Cz1.u; Cz6.y = Cz1.y;

% 7) optimal Lf
[Up2Yf_inno, Yp2Yf_inno, Uf2Yf_inno] = get_actual_matrices(A,B,C,D,K,p,f);
Lf7 = [Up2Yf_inno Yp2Yf_inno Uf2Yf_inno];
Cz7 = Lf_2_SPC(Lf7,up2usol,yp2usol,urf2usol,yrf2usol,nu,ny,p,f);
Cz7.u = Cz1.u; Cz7.y = Cz1.y;

% 8) basic IV
Lf8 = Yf_r01*Ziv8.'*pinv(Zol*Ziv8.');
% Lf8 = Yf_r01*Ziv8_2.'/(Ziv8_2*Ziv8_2.');
Cz8 = Lf_2_SPC(Lf8,up2usol,yp2usol,urf2usol,yrf2usol,nu,ny,p,f);
Cz8.u = Cz1.u; Cz8.y = Cz1.y;

% 9) approx optimal IV w/ controller information - init. basic IV
% Lf9 = Yf_r01*Ziv9.'*pinv(Zol*Ziv9.');
Lf9 = Yf_r01*Ziv9_2.'/(Ziv9_2*Ziv9_2.');
Cz9 = Lf_2_SPC(Lf9,up2usol,yp2usol,urf2usol,yrf2usol,nu,ny,p,f);
Cz9.u = Cz1.u; Cz9.y = Cz1.y;

% collect controllers and Lf estimates in cell array
nCz = 9;
[Lfs,Czs] = deal(cell(nCz,1));
Lfs{1} = Lf1; Lfs{2} = Lf2; Lfs{3} = Lf3; Lfs{4} = Lf4; Lfs{5} = Lf5;
Lfs{6} = Lf6; Lfs{7} = Lf7; Lfs{8} = Lf8; Lfs{9} = Lf9;
Czs{1} = Cz1; Czs{2} = Cz2; Czs{3} = Cz3; Czs{4} = Cz4; Czs{5} = Cz5;
Czs{6} = Cz6; Czs{7} = Cz7; Czs{8} = Cz8; Czs{9} = Cz9;

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
    Tcls{kCz} = connect(Czs{kCz},plant,Tcl_in,Tcl_out,conOpts);
end

% define data structures
[u_cl, y_cl] = deal(cell(nCz,1));
for kCz = 1:nCz
    uy_cl = lsim(Tcls{kCz},[ur1(:,f+1:end);yr1(:,f+1:end);e1],[],x1_0_cl).';
    u_cl{kCz} = uy_cl(1:nu,:);
    y_cl{kCz} = uy_cl(nu+1:end,:);
end
fprintf('Closed-loop simulations finished!\n');

%% Calculate average stage costs
fprintf('Calculate average stage costs...\n');

% average of (u_k-ur_k).' * Rk * (u_k-ur_k) over Ncl steps
cost_uk = @(u) reshape(ur1(:,1:Ncl)-u,[],1).'*kron(speye(Ncl),Rk)*reshape(ur1(:,1:Ncl)-u,[],1); 
cost_u1 = cellfun(@(u) cost_uk(u)/Ncl, u_cl);

% average of du_k.' * dRk * du_k over Ncl steps
cost_duk = @(u) reshape(u-[u0(:,end) u(:,2:end)],[],1).'*kron(speye(Ncl),dRk)*reshape(u-[u0(:,end) u(:,2:end)],[],1);
cost_u2 = cellfun(@(u) cost_duk(u)/Ncl,u_cl);

% average of (y_k-yr_k).' * Rk * (y_k-yr_k) over Ncl steps
cost_yk = @(y) reshape(yr1(:,1:Ncl)-y,[],1).'*kron(speye(Ncl),Qk)*reshape(yr1(:,1:Ncl)-y,[],1);
cost_y  = cellfun(@(y) cost_yk(y)/Ncl,y_cl);

% total costs
cost_u  = cost_u1 + cost_u2;
cost_tot= cost_u + cost_y; % total cost

%% Identification error analysis
fprintf('Calculate identification errors...\n');

FroIDerror = zeros(kCz,3);
for kCz = 1:nCz
    IDerror = Lfs{kCz}-Lf7;
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
        'u_iv4','y_iv4','u_iv5','y_iv5','u_iv9','y_iv9',...
        'cost_u1','cost_u2','cost_u','cost_y','cost_tot',...
        'FroIDerror');
    fprintf('File saved successfully!\n');
end

%% Plotting
if opts.plot
    close all;
% ===================== initial closed-loop data (& IVs) ==================
    % free (noisy) response
    [y0_3,~,~] = lsim(plant,[zeros(nu,Nbar);e0],[]); y0_3 = y0_3.';

    figure(1);
    tiledlayout(2,1,'TileSpacing','compact');
    
    ax1_11 = nexttile;
    plot(y0,'LineWidth',2); hold on; 
    plot(y0_2,'--','LineWidth',2);
    plot(y_iv4);
    plot(y_iv5);
    plot(y_iv9);
    plot(yr0);
    plot(e0);
    yLim1 = ax1_11.YLim;
    plot(y0_3); % free response
    ylim(yLim1);
    legend({'$y$','$y_{iv}^*$','$y_{iv}^{nc}$','$y_{iv}^{c,lcf}$','$y_{iv}^{c,b}$','ref','$e$','$y_{\mathrm{free}}$'},'Interpreter','latex');
    grid on;
    ylabel('$y_k$','Interpreter','latex');
    
    ax1_12 = nexttile;
    plot(u0,'LineWidth',2); hold on;
    plot(u0_2,'--','LineWidth',2);
    plot(u_iv4);
    plot(u_iv5);
    plot(u_iv9);
    grid on;
    ylabel('$u_k$','Interpreter','latex');
    xlabel('Time','Interpreter','latex');
    
    linkaxes([ax1_11,ax1_12],'x');

% ======================= closed-loop simulation data =====================
    figure(2);
    tiledlayout(2,1,'TileSpacing','compact');
    
    ax2_1 = nexttile;
    stairs(yr1, 'k-','LineWidth',2); hold on;
    for kCz = 1:nCz
        switch kCz
            case {2,7}
                plot(y_cl{kCz}, 'LineWidth',2);
            case 1
                plot(y_cl{kCz},'--');
            case {3,5}
                plot(y_cl{kCz},'--^','MarkerSize',3);
            case 4
                plot(y_cl{kCz},'--o','LineWidth',2,'MarkerSize',3);
            case 6
                plot(y_cl{kCz},'-.')
            case {8,9}
                plot(y_cl{kCz},'-','Marker','.');            
        end
    end
    ylim([-15 15]);
    leg_txt = {'$\mathcal{Z}_\mathrm{ol}$', ...
        '$\mathcal{Z}^*$', ...
        '$\mathcal{Z}_\mathrm{lcf}$',...
        '$\widehat{\mathcal{Z}}_\mathrm{nc}^*$', ...
        '$\widehat{\mathcal{Z}}_\mathrm{c,lcf}^*$', ...
        'CL-SPC', ...
        '$L_f^*$', ...
        '$\widehat{\mathcal{Z}}_\mathrm{b}$',...
        '$\widehat{\mathcal{Z}}_\mathrm{c,b}^*$',...
        'Interpreter','latex'};
    legend('ref', leg_txt{:});
    ylabel('$y_k$','Interpreter','latex')
    
    ax2_2 = nexttile;
    plot(ur1, 'k-','LineWidth',2); hold on;
    for kCz = 1:nCz
        switch kCz
            case {2,7}
                plot(u_cl{kCz}, 'LineWidth',2);
            case 1
                plot(u_cl{kCz},'--');
            case {3,5}
                plot(u_cl{kCz},'--^','MarkerSize',3);
            case 4
                plot(u_cl{kCz},'--o','LineWidth',2,'MarkerSize',3);
            case 6
                plot(u_cl{kCz},'-.')
            case {8,9}
                plot(u_cl{kCz},'-','Marker','.');            
        end
    end
    ylim([-30 30]);
    ylabel('$u_k$','Interpreter','latex');
    xlabel('Time', 'Interpreter','latex');
    
    linkaxes([ax2_1 ax2_2], 'x');
    
% ======================= visualize identification results ================
    cabsmax = 0;
    for kCz = 1:nCz
        cabsmax = max(max(abs(Lfs{kCz}),[],"all"),cabsmax);
    end
    
    figure(3);
    tiledlayout(nCz,1,'TileSpacing','compact');
    for kCz = 1:nCz
        nexttile;
        imagesc_vik(Lfs{kCz},cmax=cabsmax);
    end
    
% ======================= visualize identification error ==================
    figure(4);
    imagesc(FroIDerror);
    cmap = crameri('-davos');
    colormap(cmap);
    colorbar;

% ======================= visualize average stage costs ===================
    figure(5);
    xleg = categorical(leg_txt(1:end-2));
    xleg = reordercats(xleg,leg_txt(1:end-2));
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

function Lf = IV_GLS(Yf,Z,Phi)
    PhiZt = Phi*Z.';
    CovZPhi = (Z*Z.')\PhiZt.';
    Lf = Yf*Z.'*CovZPhi/(PhiZt*CovZPhi);
end

function Lf = IV_GLS_vec(Yf,Z,Phi,ny,fac)
    [fny,N] = size(Yf); f = fny/ny;
    PhiZt = Phi*Z.';
    YfZt = Yf*Z.';
    Lf = IV_GLS(Yf,Z,Phi); % get initial estimate
    V = Yf-Lf*Phi; % get errors
    Hf = chol(V*V.'/size(V,2),'lower');
    Re_est = Hf(1:ny,1:ny)*Hf(1:ny,1:ny).';
    Hf = Hf/kron(speye(f),Hf(1:ny,1:ny));
    Se = build_Hankel_selection_matrix_sparse(f, N);
    M1 = (Se*Se.'-speye(f*N))*fac+speye(f*N);        clear Se;
    Cov = kron(Z,Hf)*kron(M1,Re_est)*kron(Z.',Hf.'); clear M1;
    theta = kron(PhiZt,speye(fny))/Cov;              clear Cov;
    % theta = kron(PhiZt,speye(fny)) * kron(inv(Z*Z.'),Hf.'\kron(speye(f),inv(Re_est))/Hf);
    theta = (theta*kron(PhiZt.',speye(fny)))\theta*YfZt(:);
    Lf = reshape(theta,fny,[]);
end

% improve IV correlation - 2SLS
function Znew = new_IV(Z,Wp,n_exo)
    N = size(Z,2);
    if ~isempty(n_exo)
        Znew = Wp*Z.'/(Z*Z.'/N)*Z/N;
    else
        Znew = [Wp(1:n_exo,:);...
                Wp*Z.'/(Z*Z.'/N)*Z/N];
    end
end

function Se = build_Hankel_selection_matrix_sparse(f, N)
% Constructs sparse selection matrix H such that:
% vec(E_{f,N}) = H * e_tilde
% where:
% - E_{f,N} is a Hankel matrix with f rows and N columns,
% - e_tilde = [e_1; e_2; ...; e_{N+f-1}],
% - each e_t is scalar (dimension 1).

    total_rows = f * N;
    total_cols = N + f - 1;
    nz = f * N;  % number of nonzero entries (all ones)

    I = (1:nz)';                     % row indices
    J = zeros(nz, 1);                % column indices
    V = ones(nz, 1);                 % values (always 1)

    idx = 1;
    for col = 1:N
        for row = 1:f
            J(idx) = row + col - 1;
            idx = idx + 1;
        end
    end

    Se = sparse(I, J, V, total_rows, total_cols);
end

function S = newey_west(G, q)
% NEWEY_WEST Computes the Newey-West HAC covariance matrix estimate
%
% INPUTS:
%   G : m x T matrix of moment residuals (each column is one observation)
%   q : non-negative integer, bandwidth (lag truncation parameter)
%
% OUTPUT:
%   S : m x m HAC covariance matrix estimate

    [m, T] = size(G);

    % Mean-center the moment residuals
    G = G - mean(G, 2);

    % Initialize covariance matrix
    S = (G * G') / T;
    Gamma_ls = nan(m,m,q);
    % Apply weighted autocovariances
    for l = 1:q
        weight = 1 - l / (q + 1); % Bartlett kernel
        Gamma_l = (G(:, (l+1):T) * G(:, 1:(T-l))') / T;
        Gamma_ls(:,:,l) = Gamma_l;
        S = S + weight * (Gamma_l + Gamma_l');
    end
end

function Lf = IV_GMM_vec(Yf,Z,Phi,ny,nu,p,f)
    N = size(Yf,2);
    d = p*(nu+ny);
    fny = f*ny;

    % 2SLS
    Z2= new_IV(Z,Phi,d);
    Lf = Yf*Z2.'/(Z2*Z2.');
    nz = size(Z2,1);

    % get errors
    V = Yf-Lf*Phi; % get errors
    % VZ = V*Z2.';

    % intermediate values
    PhiZt = Phi*Z2.';
    YfZt  = Yf*Z2.';
    
    % estimate covariance
    Cov = newey_west(V,f);
    % Cov = kron(Z*Z.',Cov);
    % Cov = kron(speye(nz),Cov);
    theta = kron(PhiZt/(Z2*Z2.'),speye(fny)/Cov); clear Cov;
    theta = (theta*kron(PhiZt.',speye(fny)))\theta*YfZt(:);
    Lf = reshape(theta,fny,[]);
end