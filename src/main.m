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
end
[Re, p, f, N, Ncl, seed] = deal(opts.Re, opts.p, opts.f, opts.N, opts.Ncl, opts.seed);

% Requirements:
% 1) Casadi                                     v3.6.7
% 2) Control System Toolbox                     v24.2
% 3) Robust Control Toolbox                     v24.2
% 4) Statistics and Machine Learning Toolbox    v24.2
% 5) crameri_colours                            v1.09
%    Used for plotting if opts.plot = true. Obtained from
%    https://nl.mathworks.com/matlabcentral/fileexchange/68546-crameri-perceptually-uniform-scientific-colormaps

% ------------------------- add relevant paths ----------------------------
add_paths(opts);

%% Simulation settings
fprintf('Setting simulation settings...\n');

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
    cd(append('data',filesep,'raw'));
    data_dir = pwd; % -> data\raw
    cd(src_dir);

    % complete path of data directory: data\raw\<?>
    if isfield(opts,'raw_dir')
        % use data\raw\opts.raw_dir{:}
        subdir = fullfile(opts.raw_dir{:});
    else
        % use naming convention for new subdirectory
        subdir = name_subdir_data(opts);
    end
    data_dir = fullfile(data_dir,subdir);

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
    
    cCram = crameri('roma', 8); % colors

    % free (noisy) response
    [y0_free,~,~] = lsim(plant,[zeros(nu,Nbar);e0],[]); y0_free = y0_free.';
    
    figure(1);
    tiledlayout(2,1,'TileSpacing','compact');
    
    ax1_11 = nexttile; % for Yf_iv
    plot(y_iv.m6a,'Color','b','LineWidth',2); hold on;                     % reference 
    plot(y0(:,p+1:end),'LineWidth',2,'Color','k');                         % actual output
    plot(e0(:,p+1:end),'Color','r','LineStyle','--');                      % innovation noise
    plotMeanWithStd(y_iv.m2b,y_iv.s2b,Color=cCram(1,:),LineStyle='-');     % optimal IV
    plotMeanWithStd(y_iv.m4b,y_iv.s4b,Color=cCram(4,:),LineStyle='--');    % w/o controller info  
    plotMeanWithStd(u_iv.m5b,u_iv.s5b,Color=cCram(6,:),LineStyle='-');     % w/  controller info
    plotMeanWithStd(y_iv.m3a,y_iv.s3a,Color=cCram(8,:),LineStyle='-');     % LCF-IV - IV_Theta)
    yLim1 = ax1_11.YLim;
    plot(y0_free(:,p+1:end),'g');                                          % free response
    ylim(yLim1);
    legend('3a,6a) ref',     'actual output',   'noise',...
           '2b) opt. IV', '4b) w/o Cz info', '5b) w/  Cz info',...
           '3a)  LCF-IV Theta','free resp.')
    % legend({'$y$','$y_{iv}^*$','$y_{iv}^{nc}$','$y_{iv}^{c,lcf}$','$y_{iv}^{c,b}$','ref','$e$','$y_{\mathrm{free}}$'},'Interpreter','latex');
    grid on;
    ylabel('$y_k$','Interpreter','latex');
    
    ax1_12 = nexttile; % for Uf_iv
    plotMeanWithStd(u_iv.m6c,u_iv.s6c,Color='b',lineWidth=2); hold on;     % reference + 2SLS
    plot(u0(:,p+1:end),'LineWidth',2,'Color','k');                         % actual output
    plotMeanWithStd(u_iv.m2a,u_iv.s2a,Color=cCram(1,:),LineStyle='-');     % optimal IV
    plotMeanWithStd(u_iv.m2c,u_iv.s2c,Color=cCram(2,:),LineStyle='--');    % optimal IV + Yf_iv + 2SLS
    plotMeanWithStd(u_iv.m3c,u_iv.s3c,Color=cCram(3,:),LineStyle='-');     % LCF-IV + 2SLS
    plotMeanWithStd(u_iv.m4a,u_iv.s4a,Color=cCram(4,:),LineStyle='--');    % w/o controller info
    plotMeanWithStd(u_iv.m4c,u_iv.s4c,Color=cCram(5,:),LineStyle='-');     % w/o controller info + Yf_iv + 2SLS
    plotMeanWithStd(u_iv.m5a,u_iv.s5a,Color=cCram(6,:),LineStyle='-');     % w/  controller info
    plotMeanWithStd(u_iv.m5c,u_iv.s5c,Color=cCram(7,:),LineStyle='--');    % w/  controller info + Yf_iv + 2SLS
    legend('6c)  ref + 2SLS',        'actual input',    '2ab) opt. IV',...
           '2c)  opt. IV + 2SLS',    '3c)  LCF + 2SLS', '4ab) w/o Cz info',...
           '4c)  w/o Cz info + 2SLS','5ab) w/ Cz info', '5c)  w/ Cz info + 2SLS')
    grid on;
    ylabel('$u_k$','Interpreter','latex');
    xlabel('Time','Interpreter','latex');
    
    linkaxes([ax1_11,ax1_12],'x');

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
function [hFill,hLine] = plotMeanWithStd(mean_vals, std_devs, opts)
% plotMeanWithStd - Plot mean values with shaded standard deviation bands
%
% Syntax:
%   plotMeanWithStd(mean_vals, std_devs)
%   plotMeanWithStd(mean_vals, std_devs, opts)
%
% Description:
%   This function visualizes a mean curve with its uncertainty (standard
%   deviation) as a shaded region. The shaded patch is excluded from the
%   legend, ensuring the legend only describes the mean line. Plot
%   appearance can be customized via the options structure `opts`.
%
% Inputs:
%   mean_vals (1,:) double
%       Vector of mean values to plot.
%
%   std_devs (1,:) double
%       Vector of standard deviations corresponding to `mean_vals`.
%       Must be the same length as `mean_vals`.
%
%   opts (optional) struct with fields:
%       .Color  - Line/patch color. Can be:
%                       * MATLAB short color char (e.g. 'b','r','g')
%                       * RGB triplet (e.g. [0 0.447 0.741])
%                     Default: 'b'
%       .LineStyle - Line style for mean curve (e.g. '-', '--', ':')
%                     Default: '-'
%       .faceAlpha - Transparency of shaded area (0 = transparent,
%                     1 = opaque). Default: 0.2
%       .lineWidth - Width of mean line. Default: 2
%       .x         - x-axis values corresponding to `mean_vals`.
%                     Must be same length as `mean_vals`.
%                     Default: 1:length(mean_vals)
%
% Outputs:
%   hFill - handle to shaded patch (excluded from legend)
%   hLine - handle to mean line (included in legend)
%
% Example:
%   x = linspace(0,2*pi,100);
%   y = sin(x);
%   s = 0.1*ones(size(y));
%   opts.Color = [0.2 0.6 0.8];
%   opts.LineStyle = '--';
%   opts.lineWidth = 3;
%   plotMeanWithStd(y, s, opts);

    arguments
        mean_vals (1,:) double
        std_devs  (1,:) double {mustBeEqualLength(mean_vals,std_devs)}
        opts.Color  = 'b'
        opts.LineStyle = '-'
        opts.faceAlpha (1,1) double {mustBeGreaterThanOrEqual(opts.faceAlpha,0),mustBeLessThanOrEqual(opts.faceAlpha,1)} = 0.2
        opts.lineWidth (1,1) double {mustBePositive} = 2
        opts.x (1,:) double {mustBeEqualLength(opts.x,std_devs)} = 1:length(mean_vals)
    end

    x = opts.x;

    % Upper and lower bounds
    upper = mean_vals + std_devs;
    lower = mean_vals - std_devs;

    % Convert to RGB if short color char is given
    if ischar(opts.Color) || isstring(opts.Color)
        colorRGB = getColorFromChar(opts.Color);
    else
        colorRGB = opts.Color;
    end

    % Fill area between bounds (not shown in legend)
    hFill = fill([x fliplr(x)], [upper fliplr(lower)], colorRGB, ...
        'FaceAlpha', opts.faceAlpha, 'EdgeColor', 'none', ...
        'HandleVisibility','off'); 
    hold on;

    % Plot mean line (included in legend)
    hLine = plot(x, mean_vals, 'Color', colorRGB, ...
        'LineWidth', opts.lineWidth, 'LineStyle', opts.LineStyle);
end

function rgb = getColorFromChar(c)
% getColorFromChar - Convert MATLAB short color names to RGB triplets
    switch char(c)
        case 'y', rgb = [1 1 0];
        case 'm', rgb = [1 0 1];
        case 'c', rgb = [0 1 1];
        case 'r', rgb = [1 0 0];
        case 'g', rgb = [0 1 0];
        case 'b', rgb = [0 0 1];
        case 'w', rgb = [1 1 1];
        case 'k', rgb = [0 0 0];
        otherwise
            error('Unknown color specifier: %s', c);
    end
end

function mustBeEqualLength(a,b)
% mustBeEqualLength - Validation function for equal-length vectors
    if length(a) ~= length(b)
        error('Inputs must be the same length.');
    end
end

function subdir = name_subdir_data(opts)
    [Re, p, f, N, Ncl, dRk, Rk, Qk] = deal(opts.Re, opts.p, opts.f, opts.N, opts.Ncl, opts.dRk, opts.Rk, opts.Qk);

    % Helper function to trim to minimal digits in scientific notation
    trimmed_exp = @(x) regexprep(sprintf('%e', x), '(\.\d*?)0+(e[+-]?\d+)', '$1$2'); % trims trailing 0s
    % Also remove . if nothing follows
    trimmed_exp = @(x) regexprep(trimmed_exp(x), '\.(e)', '$1');

    % Apply formatting
    Re_str  = trimmed_exp(Re); N_str   = trimmed_exp(N);  Ncl_str = trimmed_exp(Ncl);
    Qk_str  = trimmed_exp(Qk); Rk_str  = trimmed_exp(Rk); dRk_str = trimmed_exp(dRk);
    subdir = sprintf('Re_%s_p_%d_f_%d_N_%s_Ncl_%s_Qk_%s_Rk_%s_dRk_%s',Re_str,p,f,N_str,Ncl_str,Qk_str,Rk_str,dRk_str);
    subdir = replace(subdir,'.','p');
    subdir = replace(subdir,'+','');
end