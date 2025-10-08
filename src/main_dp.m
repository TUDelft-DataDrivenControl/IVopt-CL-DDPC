%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           DDPC using an Optimal-IV
%           Authors: R. Dinkla, T. Oomen, J.W. van Wingerden
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opts = init_opts(N=139);
function opts = init_opts(opts)
arguments
opts.Re   (1,1) double  = 1e-1;  % innovation noise variance
opts.plot       logical = false;
opts.N    (1,1) double  = 1e3;      % number of Hankel data matrix columns
opts.f    (1,1) double  = 20;
opts.Ncl  (1,1) double  = 1500;     % simulation length of SPC
opts.dRk  (1,1) double  = 1;        % weights
opts.Rk   (1,1) double  = 1;
opts.Qk   (1,1) double  = 1e2;
opts.save       logical = true;     % save data
opts.sys  (1,1) double = 1;         % flag for model selection
end
end
[Re, N, f, Ncl] = deal(opts.Re, opts.N, opts.f, opts.Ncl);

% Requirements:
% 1) Casadi                                     v3.6.7
% 2) Control System Toolbox                     v24.2
% 3) Robust Control Toolbox                     v24.2
% 4) Statistics and Machine Learning Toolbox    v24.2
% 5) crameri_colours                            v1.09
%    Used for plotting if opts.plot = true. Obtained from
%    https://nl.mathworks.com/matlabcentral/fileexchange/68546-crameri-perceptually-uniform-scientific-colormaps
% 6) Parallel Computing Toolbox                 v24.2

% ------------------------- add relevant paths ----------------------------
add_paths(opts);

%% Simulation settings
fprintf('Setting simulation settings...\n');

% ================== initialize simulations ===============================
% get plant, initial controller (Cz0), CL-system (Tcl0), signal names, etc.
[plant,nu,ny,Cz0,Tcl0,opts,sigs] = init_sims(opts);

% ============ set p values to iterate over & number of seeds per p =======
% set p values to iterate over
pmin = max(ss2lag(plant),ss2lag(Cz0)); % take max -> if rho > p approx_IV methods deliver shorter IVs
pmax = 50;
nP   = 2;  % number of p values to iterate over
p_all = ceil(linspace(pmin,pmax,nP));

% set seeds to use for iterations
spP = 1;
seeds = reshape(1:nP*spP,spP,nP);

% ================== saving data and settings =============================
% saving this data in data\raw\sys#\dp\<subdir1>
src_dir = pwd;
cd('..'); proj_dir = pwd;
sys_dir = fullfile(pwd,'data','raw',sprintf('sys%d',opts.sys));  % -> data\raw\sys#
if ~isfolder(sys_dir)
    mkdir(sys_dir);
end
cd(src_dir);

% create data\raw\sys#\dp if it doesn't exist yet
if ~isfolder(fullfile(sys_dir,'dp'))
    mkdir(fullfile(sys_dir,'dp'))
end

% create subdir1
subdir1 = name_subdir1(pmin,pmax,nP,opts); % subdir1 name
subdir1 = fullfile(sys_dir,'dp',subdir1);  % subdir1 path
mkdir(subdir1);

% copy dependent .m files to data\raw\sys#\dp\<subdir1>\mfiles
copy_dependencies(src_dir,subdir1,'main_dp.m');

% save overall settings to data\raw\sys#\dp\<subdir1>\dp_settings.mat
save(fullfile(subdir1,'dp_settings.mat'),'pmin','pmax','nP','p_all','spP','seeds','plant','nu','ny','Cz0','Tcl0','opts','sigs');

%% ========================== iterate over p and seeds ====================
if ismember('SlurmProfile1',parallel.clusterProfiles)
    myCluster = parcluster('SlurmProfile1');
    SubmitArgsTxt = SlurmSubmitArgs('dp',30,nodes=2,tpn=48,part='compute1');
    myCluster.SubmitArguments = SubmitArgsTxt;
else
    myCluster = parcluster('local');
end
nworker = myCluster.NumWorkers; % (max.) workers per node
parpool(myCluster,nworker);

opts2 = opts;
parfor ii = 1:nP*spP
sP = struct;
sP.opts = opts2;

% ------------------------------ define p run -----------------------------
[ks,iP] = ind2sub([spP,nP],ii);
p = p_all(iP);
sP.opts.p = p;

% ----------------- initial CL-sim length & reference ---------------------
[Nbar, sP.Nbar] = deal(p + f + N -1); % sim. length of initial controller
sP.yr0  = make_reference(sP.Nbar,ny); % reference of initial controller

% ---------- references for subsequent closed-loop simulations ------------
sP.yr1 = make_reference(Ncl+f,ny); % y-ref
P0  = dcgain(plant(:,1:nu));       % DC gain
sP.ur1 = P0\sP.yr1;                % u-ref

% ----------------------- save settings for run iP ------------------------
% -> to data\raw\sys#\dp\<subdir1>\<subdir2>\<iP>_settings.mat
str_iP = iN2str(iP,nP); % zero-padded <iP> based on # of decimals for nP
subdir2 = sprintf('%s_p_%d',str_iP,p); % subdir2 name
subdir2 = fullfile(subdir1,subdir2);  % subdir2 path
if ~isfolder(subdir2)
    mkdir(subdir2); % create subdir2
end
nset_sd2 = append(str_iP,'_settings.mat');
if ~isfile(fullfile(subdir2,nset_sd2))
    save(fullfile(subdir2,nset_sd2),"-fromstruct",sP) % 'opts','Nbar','yr0','yr1','ur1'
end

% ------------------------- define seed run & noise -----------------------
sE = struct('opts',sP.opts);
seed = seeds(ii); sE.opts.seed = seed;
rng(seed);
sE.e0 = mvnrnd(zeros(ny,1),Re,Nbar).'; % innovation noise
sE.e1 = mvnrnd(zeros(ny,1),Re,Ncl).';  % innovation noise

%% run simulations
[sE.opts, sE.u0, sE.y0, sE.xcl0, sE.Z, sE.Lf, sE.Cz, sE.Tcl, sE.u_cl, sE.y_cl, sE.Cases, sE.u_iv, sE.y_iv, ...
 sE.cost_u1, sE.cost_u2, sE.cost_u, sE.cost_y, sE.cost_tot, sE.FroIDerror] ...
    = run_sims(sE.opts,sigs,plant,Cz0,Tcl0,sP.yr0,sE.e0,sP.yr1,sP.ur1,sE.e1);

%% save data
fn = sprintf('seed_%d.mat',seed);
fn = fullfile(subdir2,fn);
fn_short = strrep(fn,proj_dir,'');
fprintf('Saving data to file: \n\t%s \n',fn_short);
save(fn,"-fromstruct",sE);
%'opts','e0','e1','u0','y0','xcl0',...
% 'Z','Lf','Cz','Tcl','u_cl','y_cl','u_iv','y_iv','Cases',...
% 'cost_u1','cost_u2','cost_u','cost_y','cost_tot','FroIDerror'
fprintf('File saved successfully!\n');

end

%% Helper functions
% set name of subdir 1
function subdir1 = name_subdir1(pmin,pmax,nP,opts)
[Re, N] = deal(opts.Re, opts.N);

% Helper function to trim to minimal digits in scientific notation
trimmed_exp = @(x) regexprep(sprintf('%e', x), '(\.\d*?)0+(e[+-]?\d+)', '$1$2'); % trims trailing 0s
% Also remove . if nothing follows
trimmed_exp = @(x) regexprep(trimmed_exp(x), '\.(e)', '$1');

% Apply formatting
N_s   = trimmed_exp(N);
Re_s  = trimmed_exp(Re);
subdir1 = sprintf('p_%d_%d_%d_Re_%s_N_%s',pmin,pmax,nP,Re_s,N_s);
subdir1 = replace(subdir1,'.','p');
subdir1 = replace(subdir1,'+','');
end

% get zero-padded iN
function str_iN = iN2str(iN,nN)
nDigits = ceil(log10(nN + 1));
format = ['%0', num2str(nDigits), 'd'];
str_iN = sprintf(format,iN);
end
