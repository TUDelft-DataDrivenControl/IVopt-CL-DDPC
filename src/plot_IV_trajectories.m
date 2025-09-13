function plot_IV_trajectories(u_iv,y_iv,opts,plant,u0,y0,e0,opts2)
% plot_IV_trajectories - Compare IV trajectories with mean/std or explicit bounds
%
% Inputs:
%   u_iv, y_iv : structs with fields for mean/std OR mean/lower/upper
%   u0, y0, e0 : optional signals
%   opts.useBounds (default=false) - whether to use plotAvgWithBounds
%
% Example:
%   plot_IV_trajectories(u_iv,y_iv,u0,y0,e0, useBounds=true)

    arguments
        u_iv struct
        y_iv struct
        opts struct
        plant = [];
        u0 double = []
        y0 double = []
        e0 double = []
        opts2.useBounds logical = false
    end
    % nu = size(u_iv.m2a,1);
    [nu,p,f,N] = deal(opts.nu,opts.p,opts.f,opts.N);
    Nbar = p + f + N -1;
    
    % ============================= IV trajectories ===========================
    load("pdir.mat",'pdir');
    addpath(fullfile(pdir,'bin','external','crameri_colours'))
    cCram = crameri('roma', 8); % colors

    % Free (noisy) response only if e0 was passed
    if ~isempty(e0)
        [y0_free,~,~] = lsim(plant,[zeros(nu,Nbar);e0],[]); 
        y0_free = y0_free.';
    else
        y0_free = [];
    end

    figure(1);
    tiledlayout(2,1,'TileSpacing','compact');

    % ================= Yf_iv plot =================
    ax1_11 = nexttile; 
    plot(y_iv.m6a,'Color','b','LineWidth',2,'DisplayName','3a,6a) ref'); hold on;

    if ~isempty(y0)
        plot(y0(:,p+1:end),'LineWidth',2,'Color','k','DisplayName','actual output');
    end
    if ~isempty(e0)
        plot(e0(:,p+1:end),'Color','r','LineStyle','--','DisplayName','noise');
    end

    if opts2.useBounds
        plotAvgWithBounds(y_iv.m2b,y_iv.l2b,y_iv.u2b,Color=cCram(1,:),DisplayName='2b) opt. IV');
        plotAvgWithBounds(y_iv.m4b,y_iv.l4b,y_iv.u4b,Color=cCram(4,:),LineStyle='--',DisplayName='4b) w/o Cz info');
        plotAvgWithBounds(y_iv.m5b,y_iv.l5b,y_iv.u5b,Color=cCram(6,:),DisplayName='5b) w/ Cz info');
        plotAvgWithBounds(y_iv.m3a,y_iv.l3a,y_iv.u3a,Color=cCram(8,:),DisplayName='3a) LCF-IV Theta');
    else
        plotMeanWithStd(y_iv.m2b,y_iv.s2b,Color=cCram(1,:),DisplayName='2b) opt. IV');
        plotMeanWithStd(y_iv.m4b,y_iv.s4b,Color=cCram(4,:),LineStyle='--',DisplayName='4b) w/o Cz info');
        plotMeanWithStd(y_iv.m5b,y_iv.s5b,Color=cCram(6,:),DisplayName='5b) w/ Cz info');
        plotMeanWithStd(y_iv.m3a,y_iv.s3a,Color=cCram(8,:),DisplayName='3a) LCF-IV Theta');
    end

    yLim1 = ax1_11.YLim;
    if ~isempty(y0_free)
        plot(y0_free(:,p+1:end),'g','DisplayName','free resp.');
    end
    ylim(yLim1);
    legend show
    grid on;
    ylabel('$y_k$','Interpreter','latex');

    % ================= Uf_iv plot =================
    ax1_12 = nexttile; 
    if opts2.useBounds
        plotAvgWithBounds(u_iv.m6c,u_iv.l6c,u_iv.u6c,Color='b',lineWidth=2,DisplayName='6c) ref + 2SLS');
    else
        plotMeanWithStd(u_iv.m6c,u_iv.s6c,Color='b',lineWidth=2,DisplayName='6c) ref + 2SLS');
    end
    hold on;

    if ~isempty(u0)
        plot(u0(:,p+1:end),'LineWidth',2,'Color','k','DisplayName','actual input');
    end

    if opts2.useBounds
        plotAvgWithBounds(u_iv.m2a,u_iv.l2a,u_iv.u2a,Color=cCram(1,:),DisplayName='2a) opt. IV');
        plotAvgWithBounds(u_iv.m2c,u_iv.l2c,u_iv.u2c,Color=cCram(2,:),LineStyle='--',DisplayName='2c) opt. IV + 2SLS');
        plotAvgWithBounds(u_iv.m3c,u_iv.l3c,u_iv.u3c,Color=cCram(3,:),DisplayName='3c) LCF + 2SLS');
        plotAvgWithBounds(u_iv.m4a,u_iv.l4a,u_iv.u4a,Color=cCram(4,:),LineStyle='--',DisplayName='4a) w/o Cz info');
        plotAvgWithBounds(u_iv.m4c,u_iv.l4c,u_iv.u4c,Color=cCram(5,:),DisplayName='4c) w/o Cz info + 2SLS');
        plotAvgWithBounds(u_iv.m5a,u_iv.l5a,u_iv.u5a,Color=cCram(6,:),DisplayName='5a) w/ Cz info');
        plotAvgWithBounds(u_iv.m5c,u_iv.l5c,u_iv.u5c,Color=cCram(7,:),LineStyle='--',DisplayName='5c) w/ Cz info + 2SLS');
    else
        plotMeanWithStd(u_iv.m2a,u_iv.s2a,Color=cCram(1,:),DisplayName='2a) opt. IV');
        plotMeanWithStd(u_iv.m2c,u_iv.s2c,Color=cCram(2,:),LineStyle='--',DisplayName='2c) opt. IV + 2SLS');
        plotMeanWithStd(u_iv.m3c,u_iv.s3c,Color=cCram(3,:),DisplayName='3c) LCF + 2SLS');
        plotMeanWithStd(u_iv.m4a,u_iv.s4a,Color=cCram(4,:),LineStyle='--',DisplayName='4a) w/o Cz info');
        plotMeanWithStd(u_iv.m4c,u_iv.s4c,Color=cCram(5,:),DisplayName='4c) w/o Cz info + 2SLS');
        plotMeanWithStd(u_iv.m5a,u_iv.s5a,Color=cCram(6,:),DisplayName='5a) w/ Cz info');
        plotMeanWithStd(u_iv.m5c,u_iv.s5c,Color=cCram(7,:),LineStyle='--',DisplayName='5c) w/ Cz info + 2SLS');
    end

    legend show
    grid on;
    ylabel('$u_k$','Interpreter','latex');
    xlabel('Time','Interpreter','latex');

    linkaxes([ax1_11,ax1_12],'x');
end

%% Helper functions
function [hFill,hLine] = plotMeanWithStd(mean_vals, std_devs, opts)
arguments
    mean_vals (1,:) double
    std_devs  (1,:) double {mustBeEqualLength(mean_vals,std_devs)}
    opts.Color  = 'b'
    opts.LineStyle = '-'
    opts.faceAlpha (1,1) double {mustBeGreaterThanOrEqual(opts.faceAlpha,0),mustBeLessThanOrEqual(opts.faceAlpha,1)} = 0.2
    opts.lineWidth (1,1) double {mustBePositive} = 2
    opts.x (1,:) double {mustBeEqualLength(opts.x,std_devs)} = 1:length(mean_vals)
    opts.DisplayName string = ""   % <-- new optional input
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

% Plot mean line (included in legend, with optional DisplayName)
if opts.DisplayName ~= ""
    hLine = plot(x, mean_vals, 'Color', colorRGB, ...
        'LineWidth', opts.lineWidth, 'LineStyle', opts.LineStyle, ...
        'DisplayName', opts.DisplayName);
else
    hLine = plot(x, mean_vals, 'Color', colorRGB, ...
        'LineWidth', opts.lineWidth, 'LineStyle', opts.LineStyle);
end
end

function [hFill,hLine] = plotAvgWithBounds(average_vals, lower, upper, opts)
% plotAvgWithBounds - Plot average values with shaded region between
% provided lower and upper bounds.
%
% Syntax:
%   plotAvgWithBounds(avg, lower, upper)
%   plotAvgWithBounds(avg, lower, upper, opts)
%
% Description:
%   This function visualizes an average trajectory with its uncertainty 
%   region defined explicitly by lower and upper bounds. The shaded patch 
%   is excluded from the legend, ensuring only the mean line is labeled.
%
% Inputs:
%   average_vals (1,:) double
%       Vector of average values to plot.
%
%   lower (1,:) double
%       Lower bound trajectory (same length as average_vals).
%
%   upper (1,:) double
%       Upper bound trajectory (same length as average_vals).
%
%   opts (optional) struct with fields:
%       .Color       - Line/patch color. Char or RGB triplet. Default: 'b'
%       .LineStyle   - Line style for mean curve (default: '-')
%       .faceAlpha   - Transparency of shaded area (default: 0.2)
%       .lineWidth   - Width of mean line (default: 2)
%       .x           - x-axis values (default: 1:length(average_vals))
%       .DisplayName - Legend entry for the line (default: "")
%
% Outputs:
%   hFill - handle to shaded patch (excluded from legend)
%   hLine - handle to mean line (included in legend)

    arguments
        average_vals (1,:) double
        lower        (1,:) double {mustBeEqualLength(lower,average_vals)}
        upper        (1,:) double {mustBeEqualLength(upper,average_vals)}
        opts.Color  = 'b'
        opts.LineStyle = '-'
        opts.faceAlpha (1,1) double {mustBeGreaterThanOrEqual(opts.faceAlpha,0),...
                                     mustBeLessThanOrEqual(opts.faceAlpha,1)} = 0.2
        opts.lineWidth (1,1) double {mustBePositive} = 2
        opts.x (1,:) double {mustBeEqualLength(opts.x,average_vals)} = 1:length(average_vals)
        opts.DisplayName string = ""
    end

    x = opts.x;

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

    % Plot average line (included in legend)
    if opts.DisplayName ~= ""
        hLine = plot(x, average_vals, 'Color', colorRGB, ...
            'LineWidth', opts.lineWidth, 'LineStyle', opts.LineStyle, ...
            'DisplayName', opts.DisplayName);
    else
        hLine = plot(x, average_vals, 'Color', colorRGB, ...
            'LineWidth', opts.lineWidth, 'LineStyle', opts.LineStyle);
    end
end

function rgb = getColorFromChar(c)
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
    if length(a) ~= length(b)
        error('Inputs must be the same length.');
    end
end