function main_dN(opts)
% Perform a batch of Monte Carlo simulations sweeping over N (number of Hankel matrix columns)
% Executes spN Monte Carlo simulations for each N value in N_all, saving the data
% for subsequent analysis.
%
% Syntax:
%   main_dN()
%   main_dN(Name=Value)
%
% Important parameters to be set:
% - Nmin, Nmax: range of N values to iterate over
% - nN: number of N values to iterate over
% - spN: number of seeds (& Monte Carlo simulations) per N value
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
    opts.Nmin   (1,1) double  = 100;    % minimum N value to iterate over
    opts.Nmax   (1,1) double  = 1e4;    % maximum N value to iterate over
    opts.nN     (1,1) double  = 10;     % number of N values to iterate over
    opts.spN    (1,1) double  = 100;    % number of seeds (& Monte Carlo simulations) per N value
    opts.Re     (1,1) double  = 1e-2;   % innovation noise variance
    opts.p      (1,1) double  = 20;     % past window length
    opts.f      (1,1) double  = 20;     % future window length
    opts.Ncl    (1,1) double  = 1500;   % simulation length of SPC
    opts.dRk    (1,1) double  = 1;      % weight penalizing u_k - u_{k-1}
    opts.Rk     (1,1) double  = 1;      % weight penalizing u_k - u_{r,k}
    opts.Qk     (1,1) double  = 1e2;    % weight penalizing y_k - y_{r,k}
    opts.sys    (1,1) double = 1;       % system selection (see get_sys_info.m)
    opts.ref0   (1,:) char {mustBeMember(opts.ref0,{'make','prbs'})} = 'prbs'; % type of initial reference: 'prbs' (default) or 'make'
    opts.DMCS   (1,1) double {mustBePositive,mustBeInteger} = 1; % number of samples shift between data matrix columns (1 = Hankel, f = Page) DMCS = Data Matrix Column Shift
    opts.Cases  (1,:) cell = {'all'};   % Cases to simulate: options: 'all','iv1','CL-SPC','TrPred',etc. see CaseDefinitions.m
    opts.noSimCases cell = {};          % Cases to exclude from simulation. compatible with opts.Cases = {'all'}
end
[opts.Cases,opts.CaseDescr,opts.noSimCases] = findCases2Sim(opts.Cases,opts.noSimCases); % parse Cases
if opts.Nmin > opts.Nmax % swap if Nmin > Nmax
    [opts.Nmin, opts.Nmax] = deal(opts.Nmax, opts.Nmin);
elseif opts.Nmin == opts.Nmax
    opts.nN = 1;
end
[Re, p, f, Ncl, Nmin, Nmax, nN, spN, DMCS] = deal(opts.Re, opts.p, opts.f, opts.Ncl, opts.Nmin, opts.Nmax, opts.nN, opts.spN, opts.DMCS);

rng default;

% ------------------------- add relevant paths ----------------------------
[src_dir, ~  , ~] = fileparts(which(mfilename)); % find src directory
cd(src_dir); addpath(genpath(src_dir)); % go to, and add to path src + its subdirectories
cd('..'); addpath(fullfile('bin','casadi-v3.6.7')); cd('src'); % add path to CasADi

%% Simulation settings
fprintf('Setting simulation settings...\n');

% calculate N values to iterate over
N_all = floor(logspace(log10(Nmin),log10(Nmax),nN));

% set seeds to use for iterations
seeds = reshape(1:nN*spN,spN,nN);   % matrix with seed indices for each MC simulation

%% ================== initialize simulations ===============================
% get plant, initial controller (Cz0), CL-system (Tcl0), signal names, etc.
[plant,nu,ny,Cz0,Tcl0,opts,sigs] = init_sims(opts);

p_lb = max(ss2lag(plant),ss2lag(Cz0));
if p < p_lb
    [opts.p, p] = deal(p_lb);
    warning('p must be at least %d to accomodate the lag of the plant and initial controller. Setting p to %d.', p_lb, p_lb);
end

% Generate yr0 once with maximum length to ensure consistency across runs
Nbar_max = p + f + (Nmax-1)*DMCS; % maximum simulation length
switch opts.ref0
    case 'make'
        yr0_full = make_reference(Nbar_max,ny);
    case 'prbs'
        n_bits = ceil(log2(Nbar_max + 1));
        yr0_full = idinput(2^n_bits-1,'prbs',[0 1],[-1 1]).';
        yr0_full = repmat(yr0_full(:,1:Nbar_max),ny,1);
end

%% ================== saving data and settings =============================
% saving this data in data\sys#\ref0_<>\dN\<subdir1>
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
% create data\sys#\ref0_<>\dN if it doesn't exist yet
if ~isfolder(fullfile(ref_dir,'dN'))
    mkdir(fullfile(ref_dir,'dN'))
end
cd(src_dir);

% create subdir1 under the chosen ref folder
subdir1 = name_subdir1(Nmin,Nmax,nN,opts); % subdir1 name
subdir1 = fullfile(ref_dir,'dN',subdir1);  % subdir1 path
mkdir(subdir1);

% save overall settings to data\sys#\ref0_<>\dN\<subdir1>\dN_settings.mat
save(fullfile(subdir1,'dN_settings.mat'),'Nmin','Nmax','nN','N_all','spN','seeds','plant','nu','ny','Cz0','Tcl0','opts','sigs');

%% ========================== iterate over N and seeds ====================
load(fullfile(src_dir,'SlurmSettings.mat'),'ProfileName');
if ismember(ProfileName,parallel.clusterProfiles) % on cluster?
    % packaging input variables for use in run_X_ParCluster
    vs = struct;
    [ vs.spN, vs.nN, vs.seeds, vs.N_all, vs.p, vs.f, vs.Ncl, vs.ny, vs.nu, vs.Re, vs.plant, vs.subdir1, vs.sigs, vs.Cz0, vs.Tcl0, vs.proj_dir, vs.yr0_full, vs.ProfileName] = ...
    deal(spN,    nN,    seeds,    N_all,    p,    f,    Ncl,    ny,    nu,    Re,    plant,    subdir1,    sigs,    Cz0,    Tcl0,    proj_dir,    yr0_full,    ProfileName);

    run_X_ParCluster(opts,vs,'N',MaxTasksPerJob=50,nMins=30);

else
    fprintf('using the local profile\n');
    if isempty(gcp('nocreate'))
        myCluster = parcluster('local');
        nworker = myCluster.NumWorkers; % (max.) workers per node
        parpool(myCluster,nworker);
    end
    parfor ii = 1:nN*spN
        run_N(ii,opts,spN,nN,seeds,N_all,p,f,Ncl,ny,nu,Re,plant,subdir1,sigs,Cz0,Tcl0,proj_dir,yr0_full);
    end
end

%% Local functions
% set name of subdir 1
function subdir1 = name_subdir1(Nmin,Nmax,nN,opts)
[Re, p, f, DMCS] = deal(opts.Re, opts.p, opts.f, opts.DMCS);

% Function to trim to minimal digits in scientific notation
trimmed_exp = @(x) regexprep(sprintf('%e', x), '(\.\d*?)0+(e[+-]?\d+)', '$1$2'); % trims trailing 0s
% Also remove . if nothing follows
trimmed_exp = @(x) regexprep(trimmed_exp(x), '\.(e)', '$1');

% Apply formatting
Nmin_s = trimmed_exp(Nmin);
Nmax_s = trimmed_exp(Nmax);
Re_s  = trimmed_exp(Re);
subdir1 = sprintf('N_%s_%s_%d_Re_%s_p_%d_f_%d_DMCS_%d_%s',Nmin_s,Nmax_s,nN,Re_s,p,f,DMCS,datestr(now,'yyyymmdd_HHMM'));
subdir1 = replace(subdir1,'.','p');
subdir1 = replace(subdir1,'+','');
end
end