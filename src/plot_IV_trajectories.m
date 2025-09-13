function plot_IV_trajectories(u_iv,y_iv,opts,plant,u0,y0,e0,opts2)
% plot_IV_trajectories - Compare IV trajectories with mean/std or explicit bounds
%
% Inputs:
%   u_iv, y_iv : structs with fields for mean/std OR mean/lower/upper
%   opts        : struct with dimensions {nu,p,f,N}
%   plant       : system object (optional, for free response sim)
%   u0,y0,e0    : optional signals
%   opts2.useBounds (default=false) - use plotAvgWithBounds instead of plotMeanWithStd
%   opts2.useFill   (default=true)  - toggle shading on/off
%   opts2.showFields (default=all)  - cell array of fields to plot (e.g. {'m4a','m2b'})
%
% Example:
%   plot_IV_trajectories(u_iv,y_iv,opts,plant,u0,y0,e0,struct('useBounds',true,'useFill',false))

    arguments
        u_iv struct
        y_iv struct
        opts struct
        plant = []
        u0 double = []
        y0 double = []
        e0 double = []
        opts2.useBounds logical = false
        opts2.useFill logical = true
        opts2.showFields cell = {}   % specify fields to plot
        opts2.FontSize (1,1) double = 12;
        opts2.FaceAlpha (1,1) double = 0.2;
    end

    [nu,p,f,N] = deal(opts.nu,opts.p,opts.f,opts.N);
    Nbar = p + f + N -1;

    % ============================= IV trajectories ===========================
    load("pdir.mat",'pdir');
    addpath(fullfile(pdir,'bin','external','crameri_colours'))
    cCram = crameri('roma', 8); % colors

    % Free (noisy) response only if e0 was passed
    if ~isempty(e0) && ~isempty(plant)
        [y0_free,~,~] = lsim(plant,[zeros(nu,Nbar);e0],[]); 
        y0_free = y0_free.';
    else
        y0_free = [];
    end

    % Default fields if not specified
    if isempty(opts2.showFields)
        yFields = {'m2b','m4b','m5b','m3a'};
        uFields = {'m6c','m2a','m2c','m3c','m4a','m4c','m5a','m5c'};
    else
        % only use specified subset
        yFields = intersect(opts2.showFields, fieldnames(y_iv));
        uFields = intersect(opts2.showFields, fieldnames(u_iv));
    end

    figure(1);
    tiledlayout(2,1,'TileSpacing','compact');

    % ================= Yf_iv plot =================
    ax1_11 = nexttile; 
    plot(y_iv.m6a,'Color','k','LineWidth',2,'DisplayName','3a,6a) ref'); hold on;

    if ~isempty(y0)
        plot(y0(:,p+1:end),'LineWidth',2,'Color','b','DisplayName','actual output');
    end
    if ~isempty(e0)
        plot(e0(:,p+1:end),'Color','r','LineStyle','--','DisplayName','noise');
    end

    for k = 1:numel(yFields)
        field = yFields{k};
        switch field
            case 'm2b', label = '2b) opt. IV'; col = cCram(1,:);
            case 'm4b', label = '4b) w/o Cz info'; col = cCram(4,:);
            case 'm5b', label = '5b) w/ Cz info'; col = cCram(6,:);
            case 'm3a', label = '3a) LCF-IV Theta'; col = cCram(8,:);
            otherwise, label = field; col = [0 0 0]; % fallback
        end

        if opts2.useBounds
            plotAvgWithBounds(y_iv.(field), y_iv.(['l' field(2:end)]), y_iv.(['u' field(2:end)]), ...
                Color=col, DisplayName=label, faceAlpha=opts2.FaceAlpha*opts2.useFill);
        else
            plotMeanWithStd(y_iv.(field), y_iv.(['s' field(2:end)]), ...
                Color=col, DisplayName=label, faceAlpha=opts2.FaceAlpha*opts2.useFill);
        end
    end

    yLim1 = ax1_11.YLim;
    if ~isempty(y0_free)
        plot(y0_free(:,p+1:end),'g','DisplayName','free resp.');
    end
    ylim(yLim1);
    legend show
    grid on;
    ylabel('$y_k$','Interpreter','latex','FontSize',opts2.FontSize);

    % ================= Uf_iv plot =================
    ax1_12 = nexttile; 

    for k = 1:numel(uFields)
        field = uFields{k};
        switch field
            case 'm6c', label = '6c) ref + 2SLS'; col = 'k';
            case 'm2a', label = '2a) opt. IV'; col = cCram(1,:);
            case 'm2c', label = '2c) opt. IV + 2SLS'; col = cCram(2,:);
            case 'm3c', label = '3c) LCF + 2SLS'; col = cCram(3,:);
            case 'm4a', label = '4a) w/o Cz info'; col = cCram(4,:);
            case 'm4c', label = '4c) w/o Cz info + 2SLS'; col = cCram(5,:);
            case 'm5a', label = '5a) w/ Cz info'; col = cCram(6,:);
            case 'm5c', label = '5c) w/ Cz info + 2SLS'; col = cCram(7,:);
            otherwise, label = field; col = [0 0 0];
        end

        if opts2.useBounds
            plotAvgWithBounds(u_iv.(field), u_iv.(['l' field(2:end)]), u_iv.(['u' field(2:end)]), ...
                Color=col, DisplayName=label, faceAlpha=opts2.FaceAlpha*opts2.useFill);
        else
            plotMeanWithStd(u_iv.(field), u_iv.(['s' field(2:end)]), ...
                Color=col, DisplayName=label, faceAlpha=opts2.FaceAlpha*opts2.useFill);
        end
    end

    if ~isempty(u0)
        plot(u0(:,p+1:end),'LineWidth',2,'Color','b','DisplayName','actual input');
    end

    legend show
    grid on;
    ylabel('$u_k$','Interpreter','latex','FontSize',opts2.FontSize);
    xlabel('Time','Interpreter','latex','FontSize',opts2.FontSize);

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
if opts.faceAlpha > 0
    hFill = fill([x fliplr(x)], [upper fliplr(lower)], colorRGB, ...
        'FaceAlpha', opts.faceAlpha, 'EdgeColor', 'none', ...
        'HandleVisibility','off'); 
    hold on;
end

% Plot mean line (included in legend, with optional DisplayName)
if opts.DisplayName ~= ""
    hLine = plot(x, mean_vals, 'Color', colorRGB, ...
        'LineWidth', opts.lineWidth, 'LineStyle', opts.LineStyle, ...
        'DisplayName', opts.DisplayName);
else
    hLine = plot(x, mean_vals, 'Color', colorRGB, ...
        'LineWidth', opts.lineWidth, 'LineStyle', opts.LineStyle);
end

if opts.faceAlpha <= 0
    hold on;
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
    
    if opts.faceAlpha > 0
        % Fill area between bounds (not shown in legend)
        hFill = fill([x fliplr(x)], [upper fliplr(lower)], colorRGB, ...
            'FaceAlpha', opts.faceAlpha, 'EdgeColor', 'none', ...
            'HandleVisibility','off'); 
        hold on;
    end

    % Plot average line (included in legend)
    if opts.DisplayName ~= ""
        hLine = plot(x, average_vals, 'Color', colorRGB, ...
            'LineWidth', opts.lineWidth, 'LineStyle', opts.LineStyle, ...
            'DisplayName', opts.DisplayName);
    else
        hLine = plot(x, average_vals, 'Color', colorRGB, ...
            'LineWidth', opts.lineWidth, 'LineStyle', opts.LineStyle);
    end

    if opts.faceAlpha <= 0
        hold on;
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