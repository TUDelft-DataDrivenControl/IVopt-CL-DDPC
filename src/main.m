function main(opts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           DDPC using an Optimal-IV
%           Authors: R. Dinkla, T. Oomen, J.W. van Wingerden
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
arguments (Input)
    opts.Re   (1,1) double  = 4.81e-2;  % innovation noise variance
    opts.plot       logical = false;
    opts.p    (1,1) double  = 20;       % window lengths
    opts.f    (1,1) double  = 20;
    opts.N    (1,1) double  = 1e4;      % number of data matrix columns
    opts.Ncl  (1,1) double  = 1500;     % simulation length of SPC
    opts.dRk  (1,1) double  = 1;        % weights
    opts.Rk   (1,1) double  = 1;
    opts.Qk   (1,1) double  = 1e2;
    opts.seed (1,1) double  = 1;
    opts.save       logical = true;     % save data
    opts.raw_dir    cell;               % subdirectory of raw data directory in which to save files
    opts.sys  (1,1) double = 1;         % flag for model selection
    opts.ref0 (1,:) char {mustBeMember(opts.ref0,{'make','prbs'})} = 'prbs'; % 'make' or 'prbs'
end
[Re, p, f, N, Ncl, seed] = deal(opts.Re, opts.p, opts.f, opts.N, opts.Ncl, opts.seed);

% Requirements:
% 1) Casadi                                     v3.6.7
% 2) Control System Toolbox                     v24.2
% 3) Robust Control Toolbox                     v24.2
% 4) Statistics and Machine Learning Toolbox    v24.2
% 5) System Identification Toolbox              v24.2
% 6) crameri_colours                            v1.09
%    Used for plotting if opts.plot = true. Obtained from
%    https://nl.mathworks.com/matlabcentral/fileexchange/68546-crameri-perceptually-uniform-scientific-colormaps

rng default;

% ------------------------- add relevant paths ----------------------------
add_paths(opts);

%% Simulation settings
fprintf('Setting simulation settings...\n');

% ================== initialize simulations ===============================
% get plant, initial controller (Cz0), CL-system (Tcl0), signal names, etc.
[plant,nu,ny,Cz0,Tcl0,opts,sigs] = init_sims(opts);

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
e0 = mvnrnd(zeros(ny,1),Re,Nbar).'; % innovation noise
e1 = mvnrnd(zeros(ny,1),Re,Ncl).';  % innovation noise

%% run simulations
[opts, u0, y0, xcl0, Z, Lf, Cz, Tcl, u_cl, y_cl, Cases, u_iv, y_iv, ...
 cost_u1, cost_u2, cost_u, cost_y, cost_tot, FroIDerror] ...
    = run_sims(opts,sigs,plant,Cz0,Tcl0,yr0,e0,yr1,ur1,e1);

%% Saving data
if opts.save

    %----------------------- get destination path -------------------------
    % destination: data\raw\<?>\
    % <?> is a subdirectory defined by either [opts.raw_dir] or a naming convention
    src_dir = pwd;
    cd('..'); proj_dir = pwd;
    sys_dir = fullfile(pwd,'data','raw',sprintf('sys%d',opts.sys));  % -> data\raw\sys#
    if ~isfolder(sys_dir)
        mkdir(sys_dir);
    end
    cd(src_dir);

    % complete path of data directory: data\raw\<?>
    if isfield(opts,'raw_dir')
        % use data\raw\opts.raw_dir{:}
        subdir = fullfile(opts.raw_dir{:});
    else
        % use naming convention for new subdirectory
        subdir = name_subdir_data(opts);
    end
    data_dir = fullfile(sys_dir,subdir);

    %----------------------- make destination folder ----------------------
    % if the folder does not exist, make it and store file dependencies in it
    if ~isfolder(data_dir)
        % make data folder to store data for different seeds
        mkdir(data_dir);

        % copy dependent .m files to data_dir\mfiles
        copy_dependencies(src_dir,data_dir,'main.m');
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
        'Z','Lf','Cz','Tcl','u_cl','y_cl','u_iv','y_iv','Cases',...
        'cost_u1','cost_u2','cost_u','cost_y','cost_tot',...
        'FroIDerror');
    fprintf('File saved successfully!\n');
end

%% Plotting
if opts.plot
    close all;
    nCz = numel(Cases);
% ============================= IV trajectories ===========================
    plot_IV_trajectories(u_iv,y_iv,opts,plant,u0,y0,e0,"useBounds",false)
    
% ======================= closed-loop simulation data =====================
    % Load a Crameri colormap (e.g., 'batlow')
    colors = crameri('batlow', nCz);  % nCz colors for controllers
    
    % Define line and marker styles
    lineStyles = {'-','--',':','-.'};    % line styles
    markerStyles = {'o','s','^','d','none'};    % include 'none' to skip markers
    
    figure(2);
    tiledlayout(2,1,'TileSpacing','compact');
    
    ax2_1 = nexttile;
    stairs(yr1, 'k-','LineWidth',2,'DisplayName','ref'); hold on;
    
    for kCz = 1:nCz
        Czn = Cases{kCz}; % controller name
        ls = lineStyles{mod(kCz-1,length(lineStyles))+1};
        ms = markerStyles{mod(kCz-1,length(markerStyles))+1};
        switch Czn
            case 'CLSPC'
                stairs(y_cl.(Czn), 'Color','b','LineWidth',2,'DisplayName',Czn);
            otherwise
                stairs(y_cl.(Czn), 'Color', colors(kCz,:), 'LineStyle', ls, ...
                     'Marker', ms, 'LineWidth',1.5, 'DisplayName', Czn);
        end
    end
    
    ylim([-15 15]); grid on;
    ylabel('$y_k$','Interpreter','latex');
    legend show
    
    ax2_2 = nexttile;
    stairs(ur1, 'k-','LineWidth',2,'DisplayName','ref'); hold on;
    
    for kCz = 1:nCz
        Czn = Cases{kCz}; % controller name
        ls = lineStyles{mod(kCz-1,length(lineStyles))+1};
        ms = markerStyles{mod(kCz-1,length(markerStyles))+1};
        switch Czn
            case 'CLSPC'
                stairs(u_cl.(Czn), 'Color','b','LineWidth',2,'DisplayName',Czn);
            otherwise
                stairs(u_cl.(Czn), 'Color', colors(kCz,:), 'LineStyle', ls, ...
                     'Marker', ms, 'LineWidth',1.5, 'DisplayName', Czn);
        end
    end
    
    ylim([-30 30]); grid on;
    ylabel('$u_k$','Interpreter','latex');
    xlabel('Time', 'Interpreter','latex');
    legend show
    
    linkaxes([ax2_1 ax2_2], 'x');

    
% ======================= visualize identification results ================
    cabsmax = 0;
    for kCz = 1:nCz
        Czn = Cases{kCz}; % controller name
        cabsmax = max(max(abs(Lf.(Czn)),[],"all"),cabsmax);
    end
    
    figure(3);
    tiledlayout(nCz,1,'TileSpacing','compact');
    for kCz = 1:nCz
        Czn = Cases{kCz}; % controller name
        nexttile;
        imagesc_vik(Lf.(Czn),cmax=cabsmax);
    end
    
% ======================= visualize identification error ==================
    figure(4);
    xleg = categorical(Cases);
    xleg = reordercats(xleg,Cases);
    FroIDerror2 = zeros(nCz,3);
    for kCz = 1:nCz
        nCz = Cases{kCz};
        FroIDerror2(kCz,1) = FroIDerror.(nCz).Up;
        FroIDerror2(kCz,2) = FroIDerror.(nCz).Yp;
        FroIDerror2(kCz,3) = FroIDerror.(nCz).Uf;
    end
    subplot(1,3,1);
    bar(xleg,FroIDerror2(:,1)); grid on;
    title('Up');
    subplot(1,3,2);
    bar(xleg,FroIDerror2(:,2)); grid on;
    title('Yp');
    subplot(1,3,3);
    bar(xleg,FroIDerror2(:,3)); grid on;
    title('Uf');

% ======================= visualize average stage costs ===================
    cost_u_array = struct2array(cost_u).';
    cost_y_array = struct2array(cost_y).';
    cost_tot_array = struct2array(cost_tot).';
    figure(5);
    bar(xleg,[cost_u_array cost_y_array],'stacked');
    legend({'$\mathcal{J}_u$','$\mathcal{J}_y$'},'interpreter','latex');
    ax_bar = gca;
    ax_bar.TickLabelInterpreter = 'latex';
    max_y = max(cost_tot_array(~isoutlier(cost_tot_array)),[],'all');
    ylim(ax_bar,[0 1.05*max_y]);
end
end

%% Helper functions
function subdir = name_subdir_data(opts)
    [Re, p, N] = deal(opts.Re, opts.p, opts.N);

    % Helper function to trim to minimal digits in scientific notation
    trimmed_exp = @(x) regexprep(sprintf('%e', x), '(\.\d*?)0+(e[+-]?\d+)', '$1$2'); % trims trailing 0s
    % Also remove . if nothing follows
    trimmed_exp = @(x) regexprep(trimmed_exp(x), '\.(e)', '$1');

    % Apply formatting
    Re_str  = trimmed_exp(Re); N_str   = trimmed_exp(N);
    subdir = sprintf('Re_%s_p_%d_N_%s',Re_str,p,N_str);
    subdir = replace(subdir,'.','p');
    subdir = replace(subdir,'+','');
end