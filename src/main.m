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
    opts.N    (1,1) double  = 1e4;      % number of Hankel matrix columns
    opts.Ncl  (1,1) double  = 1500;     % simulation length of SPC
    opts.dRk  (1,1) double  = 1;        % weights
    opts.Rk   (1,1) double  = 1;
    opts.Qk   (1,1) double  = 1e2;
    opts.seed (1,1) double  = 1;
    opts.save       logical = true;     % save data
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

% go to save, and add to path raw data directory
cd(src_dir); cd('..'); cd(append('data',filesep,'raw')); raw_data_dir = pwd;
addpath(raw_data_dir);

% go back to src directory
cd(src_dir);

%% Simulation settings
fprintf('Setting simulation settings...\n');

% plant model
[plant,nx,nu,ny,A,B,C,D,K,~] = model_Landau1995();

% naming signals
uk_name = arrayfun(@(j) sprintf('u0_%d', j), 1:nu, 'UniformOutput', false);
ek_name = arrayfun(@(j) sprintf('e0_%d', j), 1:ny, 'UniformOutput', false);
yk_name = arrayfun(@(j) sprintf('y0_%d', j), 1:ny, 'UniformOutput', false);
plant.u(1:nu)     =  uk_name;
plant.u(nu+1:end) =  ek_name;
plant.y = yk_name;

% ----------------------- make/load initial controller --------------------
if isfile('Cz0_Landau1995.mat')
    Cz0 = load('Cz0_Landau1995.mat').Cz0;
else
    W1 = makeweight(33,5,0.5);  W1 = c2d(W1,plant.Ts,'tustin');
    W3 = makeweight(0.5,20,20); W3 = c2d(W3,plant.Ts,'tustin');
    [Cz0,~,~,~] = mixsyn(plant(:,1:nu),W1,[],W3);
    save('Cz0_Landau1995.mat',"Cz0");
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
yr0 = idinput(Nbar,'prbs',[],[-1 1]).';
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

% ------------------------ simulation with noise --------------------------
e0 = mvnrnd(zeros(ny,1),Re,Nbar).'; % create innovation noise

% simulate system with noise
[uy0,~,xcl0] = lsim(Tcl0,[yr0;e0],[]); uy0 = uy0.'; xcl0 = xcl0.';
u0 = uy0(1:nu,:); y0 = uy0(nu+1:end,:); clear uy0; % get inputs and outputs
xcl0_plant = xcl0(size(Cz0.A,1)+1:size(Cz0.A,1)+nx,end);  % get final state of plant

% create Hankel matrices
[~,Up_r01,Uf_r01] = make_Hankel(u0,p,f); % - data w/ noise
[~,Yp_r01,Yf_r01] = make_Hankel(y0,p,f);

% ---------------------- simulation without noise -------------------------
[uy0_2,~,~] = lsim(Tcl0,[yr0;zeros(ny,Nbar)],[]); uy0_2 = uy0_2.';
u0_2 = uy0_2(1:nu,:); y0_2 = uy0_2(nu+1:end,:); clear uy0_2;

% create Hankel matrices
[~,Up_r02,Uf_r02] = make_Hankel(u0_2,p,f); % - data w/o noise
[~,Yp_r02,~] = make_Hankel(y0_2,p,f);

%% Get Subspace Predictive Controllers
% ----------------------------- get IVs -----------------------------------
fprintf('Obtaining IVs...\n');

% 1) open-loop IV (i.e. least-squares regression)
W1 = [Up_r01;Yp_r01;Uf_r01];

% 2) optimal IV
% Zopt = [Up_r02;Yp_r02;Uf_r02];
Zopt = [Up_r01;Yp_r01;Uf_r02];

% 3) composed IV using LCF
[~,Vc,Uc] = lncf(ss(Cz0));
Hcv = make_blk_tril_toeplitz(Vc.A,Vc.B,Vc.C,Vc.D,f);
Hcu = make_blk_tril_toeplitz(Uc.A,Uc.B,Uc.C,Uc.D,f);
IV_Theta = Hcv*Uf_r01 + Hcu*Yf_r01;
Ziv3 = [Up_r01; Yp_r01; IV_Theta; Rf_yr0];

% 4) approximation of optimal IV w/o controller information
[u_iv4,y_iv4] = approx_IV_no_controller_info(u0,y0,yr0,p,p,f);
[~,Up_iv4,Uf_iv4] = make_Hankel(u_iv4,p,f);
[~,Yp_iv4,~] = make_Hankel(y_iv4,p,f);
% Ziv4 = [Up_iv4;Yp_iv4;Uf_iv4];
Ziv4 = [Up_r01;Yp_r01;Uf_iv4];

% 5) approximation of optimal IV w/ controller information
[u_iv5,y_iv5,xCz_iv] = approx_IV_controller_info(u0,y0,yr0,p,f,Ziv3,Cz0,1);
[~,Up_iv5,Uf_iv5] = make_Hankel(u_iv5,p,f);
[~,Yp_iv5,~] = make_Hankel(y_iv5,p,f);
% Ziv5 = [Up_iv5;Yp_iv5;Uf_iv5];
Ziv5 = [Up_r01;Yp_r01;Uf_iv5];

% ---------------------- get (CL-)SPC controllers -------------------------
fprintf('Obtaining controllers...\n');
[up2usol, yp2usol, urf2usol, yrf2usol] = get_solver(nu,ny,p,f,Q,R,dR);

% for initial state of SPC controllers
up0 = u0(:,end-p+1:end); up0 = up0(:);  % past input data
yp0 = y0(:,end-p+1:end); yp0 = yp0(:);  % past output data
urf = ur1(:,1:f); urf = urf(:);         % future input references
yrf = yr1(:,1:f); yrf = yrf(:);         % future output references

% 1) SPC using open-loop IV
Lf1 = Yf_r01*W1.'*pinv(W1*W1.');
Cz1 = Lf_2_SPC(Lf1,up2usol,yp2usol,urf2usol,yrf2usol,nu,ny,p,f);
urkf_name = arrayfun(@(j) sprintf('ur%d_%d',f, j), 1:nu, 'UniformOutput', false);
yrkf_name = arrayfun(@(j) sprintf('yr%d_%d',f, j), 1:ny, 'UniformOutput', false);
Cz1.u(1:ny)          = yk_name;   % y_k
Cz1.u(ny+1:ny+nu)    = urkf_name; % ur_{k+f}
Cz1.u(end-ny+(1:ny)) = yrkf_name; % yr_{k+f}
Cz1.y = uk_name;                  % u_k

% 2) optimal IV
Lf2 = Yf_r01*Zopt.'*pinv(W1*Zopt.');
[Cz2,x1_0_SPC] = Lf_2_SPC(Lf2,up2usol,yp2usol,urf2usol,yrf2usol,nu,ny,p,f,up=up0,yp=yp0,urf=urf,yrf=yrf);
Cz2.u = Cz1.u; Cz2.y = Cz1.y;

% 3) SPC using LCF
Lf3 = Yf_r01*Ziv3.'*pinv(W1*Ziv3.');
Cz3 = Lf_2_SPC(Lf3,up2usol,yp2usol,urf2usol,yrf2usol,nu,ny,p,f);
Cz3.u = Cz1.u; Cz3.y = Cz1.y;

% 4) SPC using approximation of optimal IV w/o controller information
Lf4 = Yf_r01*Ziv4.'*pinv(W1*Ziv4.');
Cz4 = Lf_2_SPC(Lf4,up2usol,yp2usol,urf2usol,yrf2usol,nu,ny,p,f);
Cz4.u = Cz1.u; Cz4.y = Cz1.y;

% 5) SPC using approximation of optimal IV w/ controller information
Lf5 = Yf_r01*Ziv5.'*pinv(W1*Ziv5.');
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

% collect controllers and Lf estimates in cell array
[Lfs,Czs] = deal(cell(7,1));
Lfs{1} = Lf1; Lfs{2} = Lf2; Lfs{3} = Lf3; Lfs{4} = Lf4; Lfs{5} = Lf5; Lfs{6} = Lf6; Lfs{7} = Lf7;
Czs{1} = Cz1; Czs{2} = Cz2; Czs{3} = Cz3; Czs{4} = Cz4; Czs{5} = Cz5; Czs{6} = Cz6; Czs{7} = Cz7;
nCz = length(Czs);

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
for kCz = 1:nCz-1           % last row is zero b/c Lf7-Lf7 = 0
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
    % Helper function to trim to minimal digits in scientific notation
    trimmed_exp = @(x) regexprep(sprintf('%e', x), '(\.\d*?)0+(e[+-]?\d+)', '$1$2'); % trims trailing 0s
    % Also remove . if nothing follows
    trimmed_exp = @(x) regexprep(trimmed_exp(x), '\.(e)', '$1');
    
    % Apply formatting
    Re_str  = trimmed_exp(Re); N_str   = trimmed_exp(N);  Ncl_str = trimmed_exp(Ncl);
    Qk_str  = trimmed_exp(Qk); Rk_str  = trimmed_exp(Rk); dRk_str = trimmed_exp(dRk);
    data_dir = sprintf('Re_%s_p_%d_f_%d_N_%s_Ncl_%s_Qk_%s_Rk_%s_dRk_%s',Re_str,p,f,N_str,Ncl_str,Qk_str,Rk_str,dRk_str);
    data_dir = replace(data_dir,'.','p');
    data_dir = replace(data_dir,'+','');
    cd(raw_data_dir);
    if ~isfolder(data_dir)
        % make data folder to store data for different seeds
        mkdir(data_dir);

        % also add .txt version of executed main.m file
        copyfile(fullfile(src_dir,append(mfilename,'.m')),...
                 fullfile(raw_data_dir,data_dir,append(mfilename,'.txt')));
    end
    fn = sprintf('seed_%d.mat',seed);
    fn = fullfile(raw_data_dir,data_dir,fn);
    fprintf('Saving data to file: \n\t%s \n',fn);
    save(fn,'Re','p','f','N','Ncl','Qk','Rk','dRk','seed','Nbar',...
        'plant','Cz0','Tcl0','yr0','yr1','ur1','e0','u0','y0','xcl0',...
        'u0_2','y0_2','Lfs','Czs','Tcls','u_cl','y_cl',...
        'cost_u1','cost_u2','cost_u','cost_y','cost_tot',...
        'FroIDerror');
    fprintf('File saved successfully!\n');
    cd(src_dir);
end

%% Plotting
if opts.plot
    close all;
% ===================== initial closed-loop data (& IVs) ==================
    figure(1);
    tiledlayout(2,1,'TileSpacing','compact');
    
    ax1_11 = nexttile;
    plot(y0); hold on; 
    plot(y0_2,'LineWidth',2);
    plot(y_iv4);
    plot(y_iv5);
    plot(yr0);
    plot(e0);
    legend({'$y$','$y_{iv}^*$','$y_{iv}^{nc}$','$y_{iv}^{c}$','ref','$e$'},'interpreter','latex');
    grid on;
    ylabel('$y_k$','Interpreter','latex');
    
    ax1_12 = nexttile;
    plot(u0); hold on;
    plot(u0_2,'LineWidth',2);
    plot(u_iv4);
    plot(u_iv5);
    grid on;
    ylabel('$u_k$','Interpreter','latex');
    xlabel('Time','Interpreter','latex');
    
    linkaxes([ax1_11,ax1_12],'x');

% ======================= closed-loop simulation data =====================
    figure(2);
    tiledlayout(2,1,'TileSpacing','compact');
    
    ax2_1 = nexttile;
    plot(yr1, 'k-','LineWidth',2); hold on;
    for kCz = 1:nCz-1
        plot(y_cl{kCz},'--');
    end
    plot(y_cl{end})
    ylim([-15 15]);
    legend('ref', ...
        '$\mathcal{Z}_\mathrm{ol}$', ...
        '$\mathcal{Z}^*$', ...
        '$\mathcal{Z}_\mathrm{lcf}$',...
        '$\widehat{\mathcal{Z}}_\mathrm{nc}^*$', ...
        '$\widehat{\mathcal{Z}}_\mathrm{c}^*$', ...
        'CL-SPC', ...
        '$L_f^*$', ...
        'Interpreter','latex');
    ylabel('$y_k$','Interpreter','latex')
    
    ax2_2 = nexttile;
    plot(ur1, 'k-','LineWidth',2); hold on;
    for kCz = 1:nCz-1
        plot(u_cl{kCz},'--');
    end
    plot(u_cl{end})
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
end

%% Helper functions
function Lf = get_Lf_CL_SPC(u1,y1,p,f,nu,ny)
    [~,Up,Uf] = make_Hankel(u1,p,1);
    [~,Yp,Yf] = make_Hankel(y1,p,1);
    L1 = Yf*pinv([Up;Yp;Uf]);
    
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
end