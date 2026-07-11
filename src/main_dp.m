function main_dp(opts)
% Perform a batch of Monte Carlo simulations sweeping over p (past window length)
% Executes spP Monte Carlo simulations for each p value in p_all, saving the data
% for subsequent analysis.
%
% Syntax:
%   main_dp()
%   main_dp(Name=Value)
%
% Important parameters to be set:
% - pmin, pmax: range of p values to iterate over
% - nP: number of p values to iterate over
% - spP: number of seeds (& Monte Carlo simulations) per p value
% - sys: determines system for which to run simulations (see get_sys_info.m, default: sys = 1)
% - other fields of opts struct (see arguments block below)
%
% Requirements:
% 1) Casadi                                     v3.6.7
% 2) Control System Toolbox                     v24.2
% 3) Robust Control Toolbox                     v24.2
% 4) Statistics and Machine Learning Toolbox    v24.2
% 5) System Identification Toolbox              v24.2
% 6) Parallel Computing Toolbox                 v24.2

arguments
    opts.pmin   (1,1) double  = 1;     % minimum p value to iterate over (auto-computed if not set)
    opts.pmax   (1,1) double  = 50;    % maximum p value to iterate over
    opts.nP     (1,1) double  = 10;    % number of p values to iterate over
    opts.spP    (1,1) double  = 100;   % number of seeds (& Monte Carlo simulations) per p value
    opts.Re     (1,1) double  = 1e-2;  % innovation noise variance
    opts.N      (1,1) double  = 1e3;   % number of Hankel data matrix columns
    opts.f      (1,1) double  = 20;    % future window length
    opts.Ncl    (1,1) double  = 1500;  % simulation length of SPC
    opts.dRk    (1,1) double  = 1;     % weight penalizing u_k - u_{k-1}
    opts.Rk     (1,1) double  = 1;     % weight penalizing u_k - u_{r,k}
    opts.Qk     (1,1) double  = 1e2;   % weight penalizing y_k - y_{r,k}
    opts.sys    (1,1) double = 1;      % system selection (see get_sys_info.m)
    opts.ref0   (1,:) char {mustBeMember(opts.ref0,{'make','prbs'})} = 'prbs'; % type of initial reference: 'prbs' (default) or 'make'
    opts.DMCS   (1,1) double {mustBePositive,mustBeInteger} = 1; % number of samples shift between data matrix columns (1 = Hankel, >1 = Page) DMCS = Data Matrix Column Shift
    opts.Cases  (1,:) cell = {'all'};   % Cases to simulate: options: 'all','iv1','CL-SPC','TrPred',etc. see CaseDefinitions.m
    opts.noSimCases cell = {};          % Cases to exclude from simulation. compatible with opts.Cases = {'all'}
end
[opts.Cases,opts.CaseDescr,opts.noSimCases] = findCases2Sim(opts.Cases,opts.noSimCases); % parse Cases
if opts.pmin > opts.pmax % swap if pmin > pmax
    [opts.pmin, opts.pmax] = deal(opts.pmax, opts.pmin);
end
[Re, N, f, Ncl, pmin, pmax, nP, spP, DMCS] = deal(opts.Re, opts.N, opts.f, opts.Ncl, opts.pmin, opts.pmax, opts.nP, opts.spP, opts.DMCS);

rng default;

% ------------------------- add relevant paths ----------------------------
[src_dir, ~  , ~] = fileparts(which(mfilename)); % find src directory
cd(src_dir); addpath(genpath(src_dir)); % go to, and add to path src + its subdirectories
cd('..'); addpath(fullfile('bin','casadi-v3.6.7')); cd('src'); % add path to CasADi

%% Simulation settings
fprintf('Setting simulation settings...\n');

% ------------------------ initialize simulations ---------------------------
% get plant, initial controller (Cz0), CL-system (Tcl0), signal names, etc.
[plant,nu,ny,Cz0,Tcl0,opts,sigs] = init_sims(opts);

% calculate p values to iterate over
p_lb = max(ss2lag(plant),ss2lag(Cz0));
if pmin < p_lb
    [opts.pmin, pmin] = deal(p_lb);
    warning('pmin must be at least %d to accomodate the lag of the plant and initial controller. Setting pmin to %d.', p_lb, p_lb);
end
if pmax < p_lb
    [opts.pmax, pmax] = deal(p_lb);
    warning('pmax must be at least %d to accomodate the lag of the plant and initial controller. Setting pmax to %d.', p_lb, p_lb);
end
if pmin == pmax
    [opts.nP, nP] = deal(1);
end
p_all = ceil(linspace(pmin,pmax,nP));

% set seeds to use for iterations
seeds = reshape(1:nP*spP,spP,nP);   % matrix with seed indices for each MC simulation

% Generate yr0 once with maximum length to ensure consistency across runs
Nbar_max = pmax + f + (N-1)*DMCS; % maximum simulation length
switch opts.ref0
    case 'make'
        yr0_full = make_reference(Nbar_max,ny);
    case 'prbs'
        n_bits = ceil(log2(Nbar_max + 1));
        yr0_full = idinput(2^n_bits-1,'prbs',[0 1],[-1 1]).';
        yr0_full = repmat(yr0_full(:,1:Nbar_max),ny,1);
end

%% ================== saving data and settings =============================
% saving this data in data\sys#\ref0_<>\dp\<subdir1>
src_dir = pwd;
cd('..'); proj_dir = pwd;
sys_dir = fullfile(pwd,'data',sprintf('sys%d',opts.sys));  % -> data\sys#
if ~isfolder(sys_dir)
    mkdir(sys_dir);
end
% create (or reuse) a reference folder under the system directory
ref_dir = fullfile(sys_dir,sprintf('ref0_%s',opts.ref0));
if ~isfolder(ref_dir)
    mkdir(ref_dir);
end
% create data\sys#\ref0_<>\dp if it doesn't exist yet
if ~isfolder(fullfile(ref_dir,'dp'))
    mkdir(fullfile(ref_dir,'dp'))
end
cd(src_dir);

% create subdir1 under the chosen ref folder
subdir1 = name_subdir1(pmin,pmax,nP,opts); % subdir1 name
subdir1 = fullfile(ref_dir,'dp',subdir1);  % subdir1 path
mkdir(subdir1);

% save overall settings to data\sys#\ref0_<>\dp\<subdir1>\dp_settings.mat
save(fullfile(subdir1,'dp_settings.mat'),'pmin','pmax','nP','p_all','spP','seeds','plant','nu','ny','Cz0','Tcl0','opts','sigs');


%% ========================== iterate over p and seeds ====================
load(fullfile(src_dir,'SlurmSettings.mat'),'ProfileName');
if ismember(ProfileName,parallel.clusterProfiles) % on cluster?
    % packaging input variables for use in run_X_ParCluster
    vs = struct;
    [ vs.spP, vs.nP, vs.seeds, vs.p_all, vs.f, vs.N, vs.Ncl, vs.ny, vs.nu, vs.Re, vs.plant, vs.subdir1, vs.sigs, vs.Cz0, vs.Tcl0, vs.proj_dir, vs.yr0_full, vs.ProfileName] = ...
    deal(spP,    nP,    seeds,    p_all,    f,    N,    Ncl,    ny,    nu,    Re,    plant,    subdir1,    sigs,    Cz0,    Tcl0,    proj_dir,    yr0_full,    ProfileName);

    run_X_ParCluster(opts,vs,'p',MaxTasksPerJob=30);

else
    fprintf('using the local profile\n');
    if isempty(gcp('nocreate'))
        myCluster = parcluster('local');
        nworker = myCluster.NumWorkers; % (max.) workers per node
        parpool(myCluster,nworker);
    end
    for iii = 1:nP*spP
        run_p(iii,opts,spP,nP,seeds,p_all,f,N,Ncl,ny,nu,Re,plant,subdir1,sigs,Cz0,Tcl0,proj_dir,yr0_full);
    end
end

%% Local functions
% set name of subdir 1
function subdir1 = name_subdir1(pmin,pmax,nP,opts)
[Re, N, f, DMCS] = deal(opts.Re, opts.N, opts.f, opts.DMCS);

% Function to trim to minimal digits in scientific notation
trimmed_exp = @(x) regexprep(sprintf('%e', x), '(\.\d*?)0+(e[+-]?\d+)', '$1$2'); % trims trailing 0s
% Also remove . if nothing follows
trimmed_exp = @(x) regexprep(trimmed_exp(x), '\.(e)', '$1');

% Apply formatting
N_s   = trimmed_exp(N);
Re_s  = trimmed_exp(Re);
subdir1 = sprintf('p_%d_%d_%d_Re_%s_N_%s_f_%d_DMCS_%d_%s',pmin,pmax,nP,Re_s,N_s,f,DMCS,datestr(now,'yyyymmdd_HHMM'));
subdir1 = replace(subdir1,'.','p');
subdir1 = replace(subdir1,'+','');
end
end