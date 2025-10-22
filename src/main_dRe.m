%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           DDPC using an Optimal-IV
%           Authors: R. Dinkla, T. Oomen, J.W. van Wingerden
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opts = init_opts(N=1e3);
[N, p, f, Ncl] = deal(opts.N, opts.p, opts.f, opts.Ncl);

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

% ============ set Re values to iterate over & number of seeds per Re =====
% set N values to iterate over
Re_min = 1e-5;
Re_max = 1;
nRe  = 10;  % number of Re values to iterate over
Re_all = logspace(log10(Re_min),log10(Re_max),nRe);

% set seeds to use for iterations
spRe = 100;
seeds = reshape(1:nRe*spRe,spRe,nRe); % spRe x nRe

% ================== initialize simulations ===============================
% get plant, initial controller (Cz0), CL-system (Tcl0), signal names, etc.
[plant,nu,ny,Cz0,Tcl0,opts,sigs] = init_sims(opts);

% ----------------- initial CL-sim length & reference ---------------------
Nbar = p + f + N -1; % sim. length of initial controller
% create initial reference (yr0).
switch opts.ref0
    case 'make'
        yr0 = make_reference(Nbar,ny); % reference of initial controller
    case 'prbs'
        n_bits = ceil(log2(Nbar + 1));
        yr0 = idinput(2^n_bits-1,'prbs',[0 1],[-10 10]).';
        yr0 = repmat(yr0(:,1:Nbar),ny,1);
end

% ---------- references for subsequent closed-loop simulations ------------
yr1 = make_reference(Ncl+f,ny); % y-ref
P0  = dcgain(plant(:,1:nu));    % DC gain
ur1 = P0\yr1;                   % u-ref

% ================== saving data and settings =============================
% saving this data in data\raw\sys#\ref0_<>\dRe\<subdir1>
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
% create data\raw\sys#\ref0_<>\dRe if it doesn't exist yet
if ~isfolder(fullfile(ref_dir,'dRe'))
    mkdir(fullfile(ref_dir,'dRe'))
end
cd(src_dir);

% create subdir1 under the chosen ref folder
subdir1 = name_subdir1(Re_min,Re_max,nRe,opts); % subdir1 name
subdir1 = fullfile(ref_dir,'dRe',subdir1);      % subdir1 path
mkdir(subdir1);

% copy dependent .m files to data\raw\sys#\ref0_<>\dRe\<subdir1>\mfiles
copy_dependencies(src_dir,subdir1,'main_dRe.m');

% save overall settings to data\raw\sys#\ref0_<>\dRe\<subdir1>\dRe_settings.mat
save(fullfile(subdir1,'dRe_settings.mat'),'Re_min','Re_max','nRe','Re_all','spRe','seeds','plant','nu','ny','Cz0','Tcl0','opts','sigs','Nbar','yr0','yr1','ur1');

%% ========================== iterate over Re & seeds =====================
if ismember('SlurmProfile1',parallel.clusterProfiles) % on cluster?
    % packaging input variables for use in run_X_ParCluster
    vs = struct;
    [ vs.spRe, vs.nRe, vs.seeds, vs.Re_all, vs.Ncl, vs.Nbar, vs.ny, vs.plant, vs.subdir1, vs.sigs, vs.Cz0, vs.Tcl0, vs.yr0, vs.yr1, vs.ur1, vs.proj_dir] = ...
    deal(spRe,    nRe,    seeds,    Re_all,    Ncl,    Nbar,    ny,    plant,    subdir1,    sigs,    Cz0,    Tcl0,    yr0,    yr1,    ur1,    proj_dir);

    run_X_ParCluster(opts,vs,'Re',MaxTasksPerJob=30);

else
    fprintf('using the local profile');
    if isempty(gcp('nocreate'))
        myCluster = parcluster('local');
        nworker = myCluster.NumWorkers; % (max.) workers per node
        parpool(myCluster,nworker);
    end
    for ii = 1:nRe*spRe
        run_Re(ii,opts,spRe,nRe,seeds,Re_all,Ncl,Nbar,ny,plant,subdir1,sigs,Cz0,Tcl0,yr0,yr1,ur1,proj_dir)
    end
end

%% Helper functions
% set name of subdir 1
function subdir1 = name_subdir1(Re_min,Re_max,nRe,opts)
[N, p, f] = deal(opts.N, opts.p, opts.f);

% Helper function to trim to minimal digits in scientific notation
trimmed_exp = @(x) regexprep(sprintf('%e', x), '(\.\d*?)0+(e[+-]?\d+)', '$1$2'); % trims trailing 0s
% Also remove . if nothing follows
trimmed_exp = @(x) regexprep(trimmed_exp(x), '\.(e)', '$1');

% Apply formatting
Re_min_s = trimmed_exp(Re_min);
Re_max_s = trimmed_exp(Re_max);
N_s  = trimmed_exp(N);
subdir1 = sprintf('Re_%s_%s_%d_p_%d_N_%s_f_%d_%s',Re_min_s,Re_max_s,nRe,p,N_s,f,datestr(now,'yyyymmdd_HHMM'));
subdir1 = replace(subdir1,'.','p');
subdir1 = replace(subdir1,'+','');
end

function opts = init_opts(opts)
arguments
opts.plot       logical = false;
opts.p    (1,1) double  = 20;       % window lengths
opts.f    (1,1) double  = 20;
opts.N    (1,1) double  = 1e4;      % number of data matrix columns
opts.Ncl  (1,1) double  = 1500;     % simulation length of SPC
opts.dRk  (1,1) double  = 1;        % weights
opts.Rk   (1,1) double  = 1;
opts.Qk   (1,1) double  = 1e2;
opts.save       logical = true;     % save data
opts.sys  (1,1) double = 1;         % flag for model selection
opts.ref0 (1,:) char {mustBeMember(opts.ref0,{'make','prbs'})} = 'prbs'; % 'make' or 'prbs'
end
end