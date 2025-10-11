%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           DDPC using an Optimal-IV
%           Authors: R. Dinkla, T. Oomen, J.W. van Wingerden
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opts = init_opts(N=1e4);
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
nP   = 2; % number of p values to iterate over
p_all = ceil(linspace(pmin,pmax,nP));

% set seeds to use for iterations
spP = 100;
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
if ismember('SlurmProfile1',parallel.clusterProfiles) % on cluster?
    % packaging input variables for use in run_X_ParCluster
    vs = struct;
    [ vs.spP, vs.nP, vs.seeds, vs.p_all, vs.f, vs.N, vs.Ncl, vs.ny, vs.nu, vs.Re, vs.plant, vs.subdir1, vs.sigs, vs.Cz0, vs.Tcl0, vs.proj_dir] = ...
    deal(spP,    nP,    seeds,    p_all,    f,    N,    Ncl,    ny,    nu,    Re,    plant,    subdir1,    sigs,    Cz0,    Tcl0,    proj_dir);

    run_X_ParCluster(opts,vs,'p');

else
    fprintf('using the local profile');
    if isempty(gcp('nocreate'))
        myCluster = parcluster('local');
        nworker = myCluster.NumWorkers; % (max.) workers per node
        parpool(myCluster,nworker);
    end
    parfor iii = 1:nP*spP
        run_p(iii,opts,spP,nP,seeds,p_all,f,N,Ncl,ny,nu,Re,plant,subdir1,sigs,Cz0,Tcl0,proj_dir)
    end
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