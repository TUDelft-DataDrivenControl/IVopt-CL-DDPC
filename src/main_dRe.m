%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           DDPC using an Optimal-IV
%           Authors: R. Dinkla, T. Oomen, J.W. van Wingerden
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opts = init_opts(N=1e3);
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
end
end
[N, p, f, Ncl] = deal(opts.N, opts.p, opts.f, opts.Ncl);

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
yr0  = make_reference(Nbar,ny); % reference of initial controller

% ---------- references for subsequent closed-loop simulations ------------
yr1 = make_reference(Ncl+f,ny); % y-ref
P0  = dcgain(plant(:,1:nu));    % DC gain
ur1 = P0\yr1;                   % u-ref

% ================== saving data and settings =============================
% saving this data in data\raw\dRe\<subdir1>
src_dir = pwd;
cd('..'); proj_dir = pwd;
cd('data'); cd('raw'); raw_dir = pwd; % -> data\raw
cd(src_dir);

% create data\raw\dRe if it doesn't exist yet
if ~isfolder(fullfile(raw_dir,'dRe'))
    mkdir(fullfile(raw_dir,'dRe'))
end

% create subdir1
subdir1 = name_subdir1(Re_min,Re_max,nRe,opts); % subdir1 name
subdir1 = fullfile(raw_dir,'dRe',subdir1);      % subdir1 path
mkdir(subdir1);

% copy dependent .m files to data\raw\dRe\<subdir1>\mfiles
copy_dependencies(src_dir,subdir1,'main_dRe.m');

% save overall settings to data\raw\dRe\<subdir1>\dRe_settings.mat
save(fullfile(subdir1,'dRe_settings.mat'),'Re_min','Re_max','nRe','Re_all','spRe','seeds','plant','nu','ny','Cz0','Tcl0','opts','sigs','Nbar','yr0','yr1','ur1');

%% ========================== iterate over Re & seeds =====================
opts2 = opts;
parfor ii = 1:nRe*spRe
opts = opts2;
s = struct;

% ------------------------------ define Re run ----------------------------
[ks,iRe] = ind2sub([spRe,nRe],ii);
Re = Re_all(iRe);
opts.Re = Re;

% ------------------------- define seed run & noise -----------------------
seed = seeds(ii); opts.seed = seed;
rng(seed);
s.e0 = mvnrnd(zeros(ny,1),Re,Nbar).'; % innovation noise
s.e1 = mvnrnd(zeros(ny,1),Re,Ncl).';  % innovation noise

% create folder to store data in
% -> to data\raw\dRe\<subdir1>\<subdir2>\seed_<seed>.mat
str_iRe = iRe2str(iRe,nRe);
subdir2 = sprintf('%s_Re_%.2e',str_iRe,Re);
subdir2 = replace(subdir2,'.','p'); % replace . with p
subdir2 = fullfile(subdir1,subdir2);
if ~isfolder(subdir2)
    mkdir(subdir2);
end

%% run simulations
[s.opts, s.u0, s.y0, s.xcl0, s.Z, s.Lf, s.Cz, s.Tcl, s.u_cl, s.y_cl, s.Cases, s.u_iv, s.y_iv, ...
 s.cost_u1, s.cost_u2, s.cost_u, s.cost_y, s.cost_tot, s.FroIDerror] ...
    = run_sims(opts,sigs,plant,Cz0,Tcl0,yr0,s.e0,yr1,ur1,s.e1);

%% save data
% ----------------------- save data for run iRe & seed ----------------
fn = sprintf('seed_%d.mat',seed);
fn = fullfile(subdir2,fn);
fn_short = strrep(fn,proj_dir,'');
fprintf('Saving data to file: \n\t%s \n',fn_short);
% save(fn,'opts','e0','e1','u0','y0','xcl0',...
%     'Z','Lf','Cz','Tcl','u_cl','y_cl','u_iv','y_iv','Cases',...
%     'cost_u1','cost_u2','cost_u','cost_y','cost_tot','FroIDerror');
save(fn,"-fromstruct",s)
fprintf('File saved successfully!\n');

end

%% Helper functions
% set name of subdir 1
function subdir1 = name_subdir1(Re_min,Re_max,nRe,opts)
[N, p, f, Ncl, dRk, Rk, Qk] = deal(opts.N, opts.p, opts.f, opts.Ncl, opts.dRk, opts.Rk, opts.Qk);

% Helper function to trim to minimal digits in scientific notation
trimmed_exp = @(x) regexprep(sprintf('%e', x), '(\.\d*?)0+(e[+-]?\d+)', '$1$2'); % trims trailing 0s
% Also remove . if nothing follows
trimmed_exp = @(x) regexprep(trimmed_exp(x), '\.(e)', '$1');

% Apply formatting
Re_min_s = trimmed_exp(Re_min);
Re_max_s = trimmed_exp(Re_max);
N_s  = trimmed_exp(N); Ncl_s = trimmed_exp(Ncl);
Qk_s  = trimmed_exp(Qk); Rk_s  = trimmed_exp(Rk); dRk_s = trimmed_exp(dRk);
subdir1 = sprintf('Re_%s_%s_%d_sys_%d_p_%d_f_%d_N_%s_Ncl_%s_Qk_%s_Rk_%s_dRk_%s',...
                   Re_min_s,Re_max_s,nRe,opts.sys,p,f,N_s,Ncl_s, Qk_s, Rk_s, dRk_s);
subdir1 = replace(subdir1,'.','p');
subdir1 = replace(subdir1,'+','');
end

% get zero-padded iRe
function str_iRe = iRe2str(iRe,nRe)
nDigits = ceil(log10(nRe + 1));
format = ['%0', num2str(nDigits), 'd'];
str_iRe = sprintf(format,iRe);
end