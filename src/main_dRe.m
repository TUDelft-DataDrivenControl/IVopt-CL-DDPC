function main_dRe(opts)
% Perform a batch of Monte Carlo simulations sweeping over Re (innovation noise variance)
% Executes spRe Monte Carlo simulations for each Re value in Re_all, saving the data
% for subsequent analysis.
%
% Syntax:
%   main_dRe()
%   main_dRe(Name=Value)
%
% Important parameters to be set:
% - Re_min, Re_max: range of Re values to iterate over
% - nRe: number of Re values to iterate over
% - spRe: number of seeds (& Monte Carlo simulations) per Re value
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
    opts.Re_min (1,1) double  = 1e-5;   % minimum Re value to iterate over
    opts.Re_max (1,1) double  = 1e-1;   % maximum Re value to iterate over
    opts.nRe    (1,1) double  = 15;     % number of Re values to iterate over
    opts.spRe   (1,1) double  = 100;    % number of seeds (& Monte Carlo simulations) per Re value
    opts.p      (1,1) double  = 20;     % past window length
    opts.f      (1,1) double  = 20;     % future window length
    opts.N      (1,1) double  = 1e3;    % number of Hankel data matrix columns
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
if opts.Re_min > opts.Re_max % swap if Re_min > Re_max
    [opts.Re_min, opts.Re_max] = deal(opts.Re_max, opts.Re_min);
elseif opts.Re_min == opts.Re_max
    opts.nRe = 1;
end
[N, p, f, Ncl, Re_min, Re_max, nRe, spRe, DMCS] = deal(opts.N, opts.p, opts.f, opts.Ncl, opts.Re_min, opts.Re_max, opts.nRe, opts.spRe, opts.DMCS);

rng default;

% ------------------------- add relevant paths ----------------------------
[src_dir, ~  , ~] = fileparts(which(mfilename)); % find src directory
cd(src_dir); addpath(genpath(src_dir)); % go to, and add to path src + its subdirectories
cd('..'); addpath(fullfile('bin','casadi-v3.6.7')); cd('src'); % add path to CasADi

%% Simulation settings
fprintf('Setting simulation settings...\n');

% calculate Re values to iterate over
Re_all = logspace(log10(Re_min),log10(Re_max),nRe);

% set seeds to use for iterations
seeds = reshape(1:nRe*spRe,spRe,nRe); % spRe x nRe

%% ================== initialize simulations ===============================
% get plant, initial controller (Cz0), CL-system (Tcl0), signal names, etc.
[plant,nu,ny,Cz0,Tcl0,opts,sigs] = init_sims(opts);

p_lb = max(ss2lag(plant),ss2lag(Cz0));
if p < p_lb
    [opts.p, p] = deal(p_lb);
    warning('p must be at least %d to accomodate the lag of the plant and initial controller. Setting p to %d.', p_lb, p_lb);
end

% ----------------- initial CL-sim length & reference ---------------------
Nbar = p + f + (N-1)*DMCS; % determine sim. length of initial controller


% create initial reference (yr0).
switch opts.ref0
    case 'make'
        yr0 = make_reference(Nbar,ny); % reference of initial controller
    case 'prbs'
        n_bits = ceil(log2(Nbar + 1));
        yr0 = idinput(2^n_bits-1,'prbs',[0 1],[-1 1]).';
        yr0 = repmat(yr0(:,1:Nbar),ny,1);
end

% ---------- references for subsequent closed-loop simulations ------------
yr1 = make_reference(Ncl+f,ny); % y-ref
P0  = dcgain(plant(:,1:nu));    % DC gain
ur1 = P0\yr1;                   % u-ref

% ================== saving data and settings =============================
% saving this data in data\sys#\ref0_<>\dRe\<subdir1>
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
% create data\sys#\ref0_<>\dRe if it doesn't exist yet
if ~isfolder(fullfile(ref_dir,'dRe'))
    mkdir(fullfile(ref_dir,'dRe'))
end
cd(src_dir);

% create subdir1 under the chosen ref folder
subdir1 = name_subdir1(Re_min,Re_max,nRe,opts); % subdir1 name
subdir1 = fullfile(ref_dir,'dRe',subdir1);      % subdir1 path
mkdir(subdir1);

% save overall settings to data\sys#\ref0_<>\dRe\<subdir1>\dRe_settings.mat
save(fullfile(subdir1,'dRe_settings.mat'),'Re_min','Re_max','nRe','Re_all','spRe','seeds','plant','nu','ny','Cz0','Tcl0','opts','sigs','Nbar','yr0','yr1','ur1');

%% ========================== iterate over Re & seeds =====================
load(fullfile(src_dir,'SlurmSettings.mat'),'ProfileName');
if ismember(ProfileName,parallel.clusterProfiles) % on cluster?
    % packaging input variables for use in run_X_ParCluster
    vs = struct;
    [ vs.spRe, vs.nRe, vs.seeds, vs.Re_all, vs.Ncl, vs.Nbar, vs.ny, vs.plant, vs.subdir1, vs.sigs, vs.Cz0, vs.Tcl0, vs.yr0, vs.yr1, vs.ur1, vs.proj_dir, vs.ProfileName] = ...
    deal(spRe,    nRe,    seeds,    Re_all,    Ncl,    Nbar,    ny,    plant,    subdir1,    sigs,    Cz0,    Tcl0,    yr0,    yr1,    ur1,    proj_dir,    ProfileName);

    run_X_ParCluster(opts,vs,'Re',MaxTasksPerJob=20);

else
    fprintf('using the local profile\n');
    if isempty(gcp('nocreate'))
        myCluster = parcluster('local');
        nworker = myCluster.NumWorkers; % (max.) workers per node
        parpool(myCluster,nworker);
    end
    parfor ii = 1:nRe*spRe
        run_Re(ii,opts,spRe,nRe,seeds,Re_all,Ncl,Nbar,ny,plant,subdir1,sigs,Cz0,Tcl0,yr0,yr1,ur1,proj_dir)
    end
end

%% Local functions
% set name of subdir 1
function subdir1 = name_subdir1(Re_min,Re_max,nRe,opts)
[N, p, f, DMCS] = deal(opts.N, opts.p, opts.f, opts.DMCS);

% Function to trim to minimal digits in scientific notation
trimmed_exp = @(x) regexprep(sprintf('%e', x), '(\.\d*?)0+(e[+-]?\d+)', '$1$2'); % trims trailing 0s
% Also remove . if nothing follows
trimmed_exp = @(x) regexprep(trimmed_exp(x), '\.(e)', '$1');

% Apply formatting
Re_min_s = trimmed_exp(Re_min);
Re_max_s = trimmed_exp(Re_max);
N_s  = trimmed_exp(N);
subdir1 = sprintf('Re_%s_%s_%d_p_%d_N_%s_f_%d_DMCS_%d_%s',Re_min_s,Re_max_s,nRe,p,N_s,f,DMCS,datestr(now,'yyyymmdd_HHMM'));
subdir1 = replace(subdir1,'.','p');
subdir1 = replace(subdir1,'+','');
end
end