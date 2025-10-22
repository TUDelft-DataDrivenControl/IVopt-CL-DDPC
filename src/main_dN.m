%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           DDPC using an Optimal-IV
%           Authors: R. Dinkla, T. Oomen, J.W. van Wingerden
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opts = init_opts(Re=1);
[Re, p, f, Ncl] = deal(opts.Re, opts.p, opts.f, opts.Ncl);

% Requirements:
% 1) Casadi                                     v3.6.7
% 2) Control System Toolbox                     v24.2
% 3) Robust Control Toolbox                     v24.2
% 4) Statistics and Machine Learning Toolbox    v24.2
% 5) System Identification Toolbox              v24.2
% 6) Parallel Computing Toolbox                 v24.2

% ------------------------- add relevant paths ----------------------------
add_paths(opts);

%% Simulation settings
fprintf('Setting simulation settings...\n');

% ============ set N values to iterate over & number of seeds per N =======
% set N values to iterate over
Nmin = 100;
Nmax = 1e4;
nN   = 10;  % number of N values to iterate over
N_all = floor(logspace(log10(Nmin),log10(Nmax),nN));

% set seeds to use for iterations
spN = 100;
seeds = reshape(1:nN*spN,spN,nN);

% ================== initialize simulations ===============================
% get plant, initial controller (Cz0), CL-system (Tcl0), signal names, etc.
[plant,nu,ny,Cz0,Tcl0,opts,sigs] = init_sims(opts);

% ================== saving data and settings =============================
% saving this data in data\raw\sys#\ref0_<>\dN\<subdir1>
src_dir = pwd;
cd('..'); proj_dir = pwd;
sys_dir = fullfile(pwd,'data','raw',sprintf('sys%d',opts.sys));  % -> data\raw\sys#
if ~isfolder(sys_dir)
    mkdir(sys_dir);
end
% create (or reuse) a reference folder under the system directory
ref_dir = fullfile(sys_dir,sprintf('ref0_%s',opts.ref0));
if ~isfolder(ref_dir)
    mkdir(ref_dir);
end
% create data\raw\sys#\ref0_<>\dN if it doesn't exist yet
if ~isfolder(fullfile(ref_dir,'dN'))
    mkdir(fullfile(ref_dir,'dN'))
end
cd(src_dir);

% create subdir1 under the chosen ref folder
subdir1 = name_subdir1(Nmin,Nmax,nN,opts); % subdir1 name
subdir1 = fullfile(ref_dir,'dN',subdir1);  % subdir1 path
mkdir(subdir1);

% copy dependent .m files to data\raw\sys#\ref0_<>\dN\<subdir1>\mfiles
copy_dependencies(src_dir,subdir1,'main_dN.m');

% save overall settings to data\raw\sys#\ref0_<>\dN\<subdir1>\dN_settings.mat
save(fullfile(subdir1,'dN_settings.mat'),'Nmin','Nmax','nN','N_all','spN','seeds','plant','nu','ny','Cz0','Tcl0','opts','sigs');

%% ========================== iterate over N and seeds ====================
if ismember('SlurmProfile1',parallel.clusterProfiles) % on cluster?
    % packaging input variables for use in run_X_ParCluster
    vs = struct;
    [ vs.spN, vs.nN, vs.seeds, vs.N_all, vs.p, vs.f, vs.Ncl, vs.ny, vs.nu, vs.Re, vs.plant, vs.subdir1, vs.sigs, vs.Cz0, vs.Tcl0, vs.proj_dir] = ...
    deal(spN,    nN,    seeds,    N_all,    p,    f,    Ncl,    ny,    nu,    Re,    plant,    subdir1,    sigs,    Cz0,    Tcl0,    proj_dir);

    run_X_ParCluster(opts,vs,'N',MaxTasksPerJob=50,nMins=30);

else
    fprintf('using the local profile');
    if isempty(gcp('nocreate'))
        myCluster = parcluster('local');
        nworker = myCluster.NumWorkers; % (max.) workers per node
        parpool(myCluster,nworker);
    end
    for ii = 1:nN*spN
        run_N(ii,opts,spN,nN,seeds,N_all,p,f,Ncl,ny,nu,Re,plant,subdir1,sigs,Cz0,Tcl0,proj_dir);
    end
end

%% Helper functions
% set name of subdir 1
function subdir1 = name_subdir1(Nmin,Nmax,nN,opts)
[Re, p, f] = deal(opts.Re, opts.p, opts.f);

% Helper function to trim to minimal digits in scientific notation
trimmed_exp = @(x) regexprep(sprintf('%e', x), '(\.\d*?)0+(e[+-]?\d+)', '$1$2'); % trims trailing 0s
% Also remove . if nothing follows
trimmed_exp = @(x) regexprep(trimmed_exp(x), '\.(e)', '$1');

% Apply formatting
Nmin_s = trimmed_exp(Nmin);
Nmax_s = trimmed_exp(Nmax);
Re_s  = trimmed_exp(Re);
subdir1 = sprintf('N_%s_%s_%d_Re_%s_p_%d_f_%d_%s',Nmin_s,Nmax_s,nN,Re_s,p,f,datestr(now,'yyyymmdd_HHMM'));
subdir1 = replace(subdir1,'.','p');
subdir1 = replace(subdir1,'+','');
end

function opts = init_opts(opts)
arguments
opts.Re   (1,1) double  = 1e-1;  % innovation noise variance
opts.plot       logical = false;
opts.p    (1,1) double  = 20;       % window lengths
opts.f    (1,1) double  = 20;
opts.Ncl  (1,1) double  = 1500;     % simulation length of SPC
opts.dRk  (1,1) double  = 1;        % weights
opts.Rk   (1,1) double  = 1;
opts.Qk   (1,1) double  = 1e2;
opts.save       logical = true;     % save data
opts.sys  (1,1) double = 1;         % flag for model selection
opts.ref0 (1,:) char {mustBeMember(opts.ref0,{'make','prbs'})} = 'prbs'; % type of reference: 'make' (default) or 'prbs'
end
end