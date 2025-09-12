%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           DDPC using an Optimal-IV
%           Authors: R. Dinkla, T. Oomen, J.W. van Wingerden
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opts = init_opts();
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
end
end
[Re, p, f, Ncl] = deal(opts.Re, opts.p, opts.f, opts.Ncl);

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
% saving this data in data\raw\dN\<subdir1>
src_dir = pwd;
cd('..'); proj_dir = pwd;
cd('data'); cd('raw'); raw_dir = pwd; % -> data\raw
cd(src_dir);

% create data\raw\dN if it doesn't exist yet
if ~isfolder(fullfile(raw_dir,'dN'))
    mkdir(fullfile(raw_dir,'dN'))
end

% create subdir1
subdir1 = name_subdir1(Nmin,Nmax,nN,opts); % subdir1 name
subdir1 = fullfile(raw_dir,'dN',subdir1);  % subdir1 path
mkdir(subdir1);

% copy dependent .m files to data\raw\dN\<subdir1>\mfiles
copy_dependencies(src_dir,subdir1,'main_dN.m');

% save overall settings to data\raw\dN\<subdir1>\dN_settings.mat
save(fullfile(subdir1,'dN_settings.mat'),'Nmin','Nmax','nN','N_all','spN','seeds','plant','nu','ny','Cz0','Tcl0','opts','sigs');

%% ========================== iterate over N and seeds ====================
opts2 = opts;
parfor ii = 1:nN*spN
sN = struct;
sN.opts = opts2;

% ------------------------------ define N run -----------------------------
[ks,iN] = ind2sub([spN,nN],ii);
N = N_all(iN);
sN.opts.N = N;

% ----------------- initial CL-sim length & reference ---------------------
sN.Nbar = p + f + N -1; % sim. length of initial controller
sN.yr0  = make_reference(sN.Nbar,ny); % reference of initial controller

% ---------- references for subsequent closed-loop simulations ------------
sN.yr1 = make_reference(Ncl+f,ny); % y-ref
P0  = dcgain(plant(:,1:nu));       % DC gain
sN.ur1 = P0\sN.yr1;                % u-ref

% ----------------------- save settings for run iN ------------------------
% -> to data\raw\dN\<subdir1>\<subdir2>\<iN>_settings.mat
str_iN = iN2str(iN,nN); % zero-padded <iN> based on # of decimals for nN
subdir2 = sprintf('%s_N_%d',str_iN,N); % subdir2 name
subdir2 = fullfile(subdir1,subdir2);  % subdir2 path
if ~isfolder(subdir2)
    mkdir(subdir2); % create subdir2
end
nset_sd2 = append(str_iN,'_settings.mat');
if ~isfile(fullfile(subdir2,nset_sd2))
    save(fullfile(subdir2,nset_sd2),"-fromstruct",sN) % 'opts','Nbar','yr0','yr1','ur1'
end

% ------------------------- define seed run & noise -----------------------
sE = struct('opts',sN.opts);
seed = seeds(ii); sE.opts.seed = seed;
rng(seed);
sE.e0 = mvnrnd(zeros(ny,1),Re,Nbar).'; % innovation noise
sE.e1 = mvnrnd(zeros(ny,1),Re,Ncl).';  % innovation noise

%% run simulations
[sE.opts, sE.u0, sE.y0, sE.xcl0, sE.Z, sE.Lf, sE.Cz, sE.Tcl, sE.u_cl, sE.y_cl, sE.Cases, sE.u_iv, sE.y_iv, ...
 sE.cost_u1, sE.cost_u2, sE.cost_u, sE.cost_y, sE.cost_tot, sE.FroIDerror] ...
    = run_sims(sE.opts,sigs,plant,Cz0,Tcl0,sN.yr0,sE.e0,sN.yr1,sN.ur1,sE.e1);

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
function subdir1 = name_subdir1(Nmin,Nmax,nN,opts)
[Re, p, f, Ncl, dRk, Rk, Qk] = deal(opts.Re, opts.p, opts.f, opts.Ncl, opts.dRk, opts.Rk, opts.Qk);

% Helper function to trim to minimal digits in scientific notation
trimmed_exp = @(x) regexprep(sprintf('%e', x), '(\.\d*?)0+(e[+-]?\d+)', '$1$2'); % trims trailing 0s
% Also remove . if nothing follows
trimmed_exp = @(x) regexprep(trimmed_exp(x), '\.(e)', '$1');

% Apply formatting
Nmin_s = trimmed_exp(Nmin);
Nmax_s = trimmed_exp(Nmax);
Re_s  = trimmed_exp(Re); Ncl_s = trimmed_exp(Ncl);
Qk_s  = trimmed_exp(Qk); Rk_s  = trimmed_exp(Rk); dRk_s = trimmed_exp(dRk);
subdir1 = sprintf('N_%s_%s_%d_sys_%d_Re_%s_p_%d_f_%d_Ncl_%s_Qk_%s_Rk_%s_dRk_%s',...
                   Nmin_s,Nmax_s,nN,opts.sys,Re_s,p,f, Ncl_s, Qk_s, Rk_s, dRk_s);
subdir1 = replace(subdir1,'.','p');
subdir1 = replace(subdir1,'+','');
end

% get zero-padded iN
function str_iN = iN2str(iN,nN)
nDigits = ceil(log10(nN + 1));
format = ['%0', num2str(nDigits), 'd'];
str_iN = sprintf(format,iN);
end
