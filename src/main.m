function main(opts)
%% Single Monte Carlo simulation wrapper with configurable settings and seed
% Executes a complete end-to-end DDPC session with specified parameters (N, p, Re, etc.)
% and random seed, returning performance metrics for instrumental variable comparison.
%
% Requirements:
% 1) Casadi                                     v3.6.7
% 2) Control System Toolbox                     v24.2
% 3) Robust Control Toolbox                     v24.2
% 4) Statistics and Machine Learning Toolbox    v24.2
% 5) System Identification Toolbox              v24.2

arguments (Input)
    opts.Re   (1,1) double  = 1e-2;  % innovation noise variance
    opts.p    (1,1) double  = 20;       % past window length
    opts.f    (1,1) double  = 20;       % future window length
    opts.N    (1,1) double  = 1e3;      % number of Hankel data matrix columns
    opts.Ncl  (1,1) double  = 1500;     % simulation length of SPC
    opts.dRk  (1,1) double  = 1;        % weight penalizing u_k - u_{k-1}
    opts.Rk   (1,1) double  = 1;        % weight penalizing u_k - u_{r,k}
    opts.Qk   (1,1) double  = 1e2;      % weight penalizing y_k - y_{r,k}
    opts.seed (1,1) double  = 1;        % random seed for reproducibility
    opts.save       logical = true;     % save data
    opts.raw_dir    cell;               % subdirectory of raw data directory in which to save files
    opts.sys    (1,1) double = 1;       % flag for model selection
    opts.ref0   (1,:) char {mustBeMember(opts.ref0,{'make','prbs'})} = 'prbs'; % 'make' or 'prbs'
    opts.DMCS   (1,1) double {mustBePositive,mustBeInteger} = 1; % number of samples shift between data matrix columns (1 = Hankel, >1 = Page) DMCS = Data Matrix Column Shift
    opts.Cases  (1,:) cell = {'all'};   % Cases to simulate: options: 'all','iv1','CL-SPC','TrPred',etc. see CaseDefinitions.m
    opts.noSimCases cell = {};          % Cases to exclude from simulation. compatible with opts.Cases = {'all'}
end
[opts.Cases,opts.CaseDescr,opts.noSimCases] = findCases2Sim(opts.Cases,opts.noSimCases); % parse Cases
[Re, p, f, N, Ncl, seed] = deal(opts.Re, opts.p, opts.f, opts.N, opts.Ncl, opts.seed);

rng default;

% ------------------------- add relevant paths ----------------------------
[src_dir, ~  , ~] = fileparts(which(mfilename)); % find src directory
cd(src_dir); addpath(genpath(src_dir)); % go to, and add to path src + its subdirectories
cd('..'); addpath(fullfile('bin','casadi-v3.6.7')); cd('src'); % add path to CasADi

%% Simulation settings
fprintf('Setting simulation settings...\n');

% ================== initialize simulations ===============================
% get plant, initial controller (Cz0), CL-system (Tcl0), signal names, etc.
[plant,nu,ny,Cz0,Tcl0,opts,sigs] = init_sims(opts);

p_lb = max(ss2lag(plant),ss2lag(Cz0));
if p < p_lb
    [opts.p, p] = deal(p_lb);
    warning('p must be at least %d to accomodate the lag of the plant and initial controller. Setting p to %d.', p_lb, p_lb);
end

% ----------------- initial CL-sim length & reference ---------------------
Nbar = p + f + N -1; % sim. length of initial controller
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

% ==================== seed-dependent from here onwards ===================
rng(seed);
e0 = mvnrnd(zeros(ny,1),Re,Nbar).'; % innovation noise for initial closed-loop simulation
e1 = mvnrnd(zeros(ny,1),Re,Ncl).';  % subsequent innovation noise

%% run simulations
[opts, u0, y0, xcl0, Z, Lf, Cz, Tcl, u_cl, y_cl] ...
    = main_MC(opts,sigs,plant,Cz0,Tcl0,yr0,e0,yr1,ur1,e1);

%% Saving data
if opts.save

    %----------------------- get destination path -------------------------
    % destination: data\sys#\ref0_<type>\<params>\
    % subdirectory structure includes ref0 type and parameter naming convention
    src_dir = pwd;
    cd('..'); proj_dir = pwd;
    sys_dir = fullfile(pwd,'data',sprintf('sys%d',opts.sys));  % -> data\sys#
    if ~isfolder(sys_dir)
        mkdir(sys_dir);
    end
    cd(src_dir);

    % complete path of data directory: data\sys#\ref0_<type>\<params>\
    if isfield(opts,'raw_dir')
        % use data\opts.raw_dir{:}
        subdir = fullfile(opts.raw_dir{:});
        data_dir = fullfile(sys_dir,subdir);
    else
        % use naming convention: data\sys#\ref0_<type>\<params>\
        ref0_dir = sprintf('ref0_%s',opts.ref0);  % ref0_make or ref0_prbs
        param_subdir = name_subdir_data(opts);
        data_dir = fullfile(sys_dir,ref0_dir,param_subdir);
    end

    %----------------------- make destination folder ----------------------
    % if the folder does not exist, make it and store file dependencies in it
    if ~isfolder(data_dir)
        % make data folder to store data for different seeds
        mkdir(data_dir);
    end
    
    %----------------------- add path to data directory -------------------
    addpath(data_dir);
    
    %----------------------- save data to file in destination -------------
    fn = sprintf('seed_%d.mat',seed);
    fn = fullfile(data_dir,fn);
    fn_short = strrep(fn,proj_dir,'');
    fprintf('Saving data to file: \n\t%s \n',fn_short);
    save(fn,'opts',...
        'plant','Cz0','Tcl0','yr0','yr1','ur1','e1','e0','u0','y0','xcl0',...
        'Z','Lf','Cz','Tcl','u_cl','y_cl');
    fprintf('File saved successfully!\n');
end
end

%% Local functions
function subdir = name_subdir_data(opts)
    [Re, N, p, f, Ncl, Qk, Rk, dRk] = deal(opts.Re, opts.N, opts.p, opts.f, opts.Ncl, opts.Qk, opts.Rk, opts.dRk);

    % Function to trim to minimal digits in scientific notation
    trimmed_exp = @(x) regexprep(sprintf('%e', x), '(\.\d*?)0+(e[+-]?\d+)', '$1$2'); % trims trailing 0s
    % Also remove . if nothing follows
    trimmed_exp = @(x) regexprep(trimmed_exp(x), '\.(e)', '$1');

    % Apply formatting for scientific notation parameters
    Re_str = trimmed_exp(Re);
    N_str = trimmed_exp(N);
    
    % Construct subdirectory name
    subdir = sprintf('Re_%s_N_%s_p_%d_f_%d_Ncl_%d_Qk_%g_Rk_%g_dRk_%g', ...
        Re_str, N_str, p, f, Ncl, Qk, Rk, dRk);
    subdir = replace(subdir,'.','p');
    subdir = replace(subdir,'+','');
end