function [ax1, fig] = Fig_IV_approx(data_type,fig_dir,Uf_ivs)
%FIG_IV_APPROX Plot log-log Frobenius norm errors for selected Uf IV approximations.
% Load processed metrics (m1.Uf*) and configuration values for the selected sweep.
fig_file = fullfile(fig_dir,'processed_data.mat');
load(fig_file);
load(fullfile(fig_dir,sprintf('d%s_settings.mat',data_type)));
[p,f,nu,ny,N] = deal(opts.p,opts.f,opts.nu,opts.ny,opts.N);

FS_Tick   = 9;
FS_Label  = 12;
Fs_Legend = 10;

FigPos = [50 50 252 325/2];
Units = 'points';
CrameriColor = 'batlow'; %'romaO'
LineWidth = 2;
MarkerSize = 6;

% Determine x-axis
switch data_type
    case 'N'
        xlab = '$N$';
    case {'Re','p'}
        if strcmp(data_type,'Re')
            xlab = '$\Sigma_e$';
            X_all = Re_all;
        else
            xlab = '$p$';
            X_all = p_all;
        end        
end

Colors = make_color_struct(Uf_ivs, CrameriColor);

% Small wrappers keep axis-label creation centralized and consistent.
set_xlabel   = @() xlabel(xlab,'interpreter','latex','FontSize',FS_Label);
set_ylabel = @(UorY) ylabel(['$\big\| \widetilde{',UorY,'}_{\mathrm{f}}-\alpha\big\|_\mathrm{F}$'],...
    'Interpreter','latex','FontSize',FS_Label*1.25);%,'Rotation',0

% Build a single axes figure for the log-log comparison plot.
fig = figure('Units',Units,'Position',FigPos);
fig.Units = 'pixels';
ax1 = axes;

% Keep only requested IV cases that are present in m1.Uf.
Uf_ivs = intersect(fieldnames(m1.Uf).',Uf_ivs);
Uf_av_fs = replace(Uf_ivs,'iv','m');
num_Uf_ivs = numel(Uf_ivs);

% Plot each selected IV variant and collect legend handles/labels.
u_iv = struct;
u_plot_handles = [];
u_labels = {};
for kIV = 1:num_Uf_ivs
    av_name = Uf_av_fs{kIV};
    iv_name = Uf_ivs{kIV};
    
    u_iv.(av_name) = m1.Uf.(iv_name).median;

    % Assign display label and style per IV variant.
    switch iv_name
        case 'iv4a', label = '$\hat{\widetilde{U}}_{\mathrm{f,4}}$'; LineStyle = '-';  Marker = 'none';
        case 'iv5a', label = '$\hat{\widetilde{U}}_{\mathrm{f,5}}$'; LineStyle = '--';  Marker = 'none';
        case 'iv4b', label = '$\hat{\widetilde{U}}_{\mathrm{f,4b}}$'; LineStyle = '--'; Marker = 'none';
        case 'iv5b', label = '$\hat{\widetilde{U}}_{\mathrm{f,5b}}$'; LineStyle = '--'; Marker = 's';
        case 'iv4c', label = '$\hat{\widetilde{U}}_{\mathrm{f,4c}}$'; LineStyle = ':';  Marker = 'o';
        case 'iv5c', label = '$\hat{\widetilde{U}}_{\mathrm{f,5c}}$'; LineStyle = ':';  Marker = 'd';
        case 'iv4d', label = '$\hat{\widetilde{U}}_{\mathrm{f,4d}}$'; LineStyle = '-.'; Marker = '^';
        case 'iv5d', label = '$\hat{\widetilde{U}}_{\mathrm{f,5d}}$'; LineStyle = '-.'; Marker = 'v';
        case 'iv1',  label = '$U_{\mathrm{f}}\vphantom{\hat{\widetilde{U}}}$'; LineStyle = '-';  Marker = 'none';
        case 'iv6',  label = '$W_{\mathrm{f}}$';                      LineStyle = '-';  Marker = 'x';
        otherwise,   label = iv_name;                                 LineStyle = '-';  Marker = 'none';
    end
    col = Colors.(iv_name);

    h = loglog(X_all, u_iv.(av_name), 'Marker', Marker, 'MarkerSize', MarkerSize,'Color', col, 'LineStyle', LineStyle, 'LineWidth', LineWidth);
    u_plot_handles = [u_plot_handles, h];
    u_labels{kIV} = label;
    hold on;
end

% Final axes styling.
set_ylabel('U');
set_xlabel();
grid on;

ax1.FontSize = FS_Tick;
ax1.YLabel.FontSize = FS_Label;
ax1.XLabel.FontSize = FS_Label;

leg = legend(ax1, u_plot_handles, u_labels, 'Interpreter', 'latex', 'FontSize', Fs_Legend,...
    'NumColumns', 2, 'Box', 'on', 'IconColumnWidth', 20,'Location','southeast');

% Add a small "$\alpha=$" title to the left of the legend entries.
vertices = [leg.Position(1:2); leg.Position(1), leg.Position(2) + leg.Position(4)];
h_title = make_alpha_title_left(Fs_Legend + 1, ax1, vertices, 'figure');

end

%% Local functions
function color_struct = make_color_struct(IV_names,CrameriColor)
%MAKE_COLOR_STRUCT Map IV names to consistent color groups.
% IV4* share one color, IV5* share one color, and iv6 is forced to black.
    nIVs = numel(IV_names);
    color_types = nan(1,nIVs);
    for k = 1:nIVs
        IV_name = IV_names{k};
        switch IV_name(3:end)
            case '1'
                col = 1;
            case '2'
                col = 2;
            case {'4a','4b','4c','4d'}
                col = 3;
            case {'5a','5b','5c','5d'}
                col = 4;
            case '6'
                col = 5;
            otherwise
                continue
        end
        color_types(k) = col;
    end
    color_types2 = unique(color_types(~isnan(color_types)));
    nColTypes = numel(color_types2);

    color_struct = struct;
    if contains(IV_names,'iv6')
        cCram = [crameri(CrameriColor, nColTypes-1); zeros(1,3)]; % colors
    else
        cCram = crameri(CrameriColor, nColTypes); % colors
    end
    for k = 1:nIVs
        if isnan(color_types(k))
            continue
        end
        color_type = color_types(k);
        idx = find(color_types2 == color_type);
        color_struct.(IV_names{k}) = cCram(idx,:);
    end
end

% position relative to axis -> position relative to figure
function RelPosFig = RelPosAx2Fig(ax, RelPosAx, varargin)
    narginchk(2,3);

    % save old units
    axUnits = ax.Units;
    fig = gcf; FigUnits = fig.Units;
    
    % set units to normalized
    ax.Units = 'normalized';
    fig.Units = 'normalized';
    
    % By default convert [x y w h]; optionally convert selected coordinate types.
    if numel(RelPosAx) == 4
        coordTypes = {'x','y','w','h'};
    elseif numel(RelPosAx) ~= numel(varargin{1})
        error('Please specify coordinate types for each element of RelPosAx (e.g. ''x'',''y'',''w'',''h'')');
    else
        coordTypes = varargin{1};
    end

    axPos = ax.Position;
    [RelPosFig.x, RelPosFig.y, RelPosFig.w, RelPosFig.h] = deal(nan);
    for k = 1:numel(coordTypes)
        switch coordTypes{k}
            case 'x'
                RelPosFig.x = axPos(1)+RelPosAx(1)*axPos(3);
            case 'y'
                RelPosFig.y = axPos(2)+RelPosAx(2)*axPos(4);
            case 'w'
                RelPosFig.w = RelPosAx(3)*axPos(3);
            case 'h'
                RelPosFig.h = RelPosAx(4)*axPos(4);
        end
    end
    RelPosFig = [RelPosFig.x, RelPosFig.y, RelPosFig.w, RelPosFig.h];
    RelPosFig = RelPosFig(~isnan(RelPosFig));

    % reset units
    fig.Units = FigUnits;
    ax.Units = axUnits;
end

% Convert position from figure frame of reference to axes frame of reference (inverse of RelPosAx2Fig)
function RelPosAx = RelPosFig2Ax(ax, RelPosFig, varargin)
    narginchk(2,3);
    
    % save old units
    axUnits = ax.Units;
    fig = gcf; FigUnits = fig.Units;
    
    % set units to normalized
    ax.Units = 'normalized';
    fig.Units = 'normalized';
    
    % By default convert [x y w h]; optionally convert selected coordinate types.
    if numel(RelPosFig) == 4
        coordTypes = {'x','y','w','h'};
    elseif numel(RelPosFig) ~= numel(varargin{1})
        error('Please specify coordinate types for each element of RelPosFig (e.g. ''x'',''y'',''w'',''h'')');
    else
        coordTypes = varargin{1};
    end

    axPos = ax.Position;
    [RelPosAx_x, RelPosAx_y, RelPosAx_w, RelPosAx_h] = deal(nan);
    for k = 1:numel(coordTypes)
        switch coordTypes{k}
            case 'x'
                RelPosAx_x = (RelPosFig(1) - axPos(1)) / axPos(3);
            case 'y'
                RelPosAx_y = (RelPosFig(2) - axPos(2)) / axPos(4);
            case 'w'
                RelPosAx_w = RelPosFig(3) / axPos(3);
            case 'h'
                RelPosAx_h = RelPosFig(4) / axPos(4);
        end
    end
    RelPosAx = [RelPosAx_x, RelPosAx_y, RelPosAx_w, RelPosAx_h];
    RelPosAx = RelPosAx(~isnan(RelPosAx));

    % reset units
    fig.Units = FigUnits;
    ax.Units = axUnits;
end

% Create text object on the first overlay axis with background box
function h_title = make_alpha_title_left(FontSize,ax,vertices,FigAxMode)
%MAKE_ALPHA_TITLE_LEFT Place "$\alpha=$" left of legend
        if strcmp(FigAxMode,'axes')
            vertices(1,:) = RelPosAx2Fig(ax,vertices(1,:),{'x','y'});
            vertices(2,:) = RelPosAx2Fig(ax,vertices(2,:),{'x','y'});
        elseif ~strcmp(FigAxMode,'figure')
            error('Either specify "axes" or "figure" as relative coordinate type.')
        end
        RelPosTitle = mean(vertices) - [0.02 0]; % Fig
        RelPosTitle  = RelPosFig2Ax(ax, RelPosTitle, {'x','y'}); % -> Ax
        h_title = text(ax, RelPosTitle(1), RelPosTitle(2), '$\alpha=$', ...
            'Units', 'normalized', ...
            'Interpreter', 'latex', ...
            'FontSize', FontSize, ...
            'HorizontalAlignment', 'right', ...
            'VerticalAlignment', 'middle', ...
            'BackgroundColor', 'white', ...
            'EdgeColor', 'none');
end