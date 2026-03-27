function [ax, fig] = Fig_IV_approx(data_type,fig_dir,legPosOnAx,legendEntriesPerColumn)
fig_file = fullfile(fig_dir,'processed_data.mat');
load(fig_file);
load(fullfile(fig_dir,sprintf('d%s_settings.mat',data_type)));
[p,f,nu,ny,N] = deal(opts.p,opts.f,opts.nu,opts.ny,opts.N);

FS_Tick   = 9;
FS_Label  = 12;
Fs_Legend = 10;

% Extract and set up optional arguments
FigPos = [50 50 252 325];
Units = 'points';
CrameriColor = 'romaO'; % batlow
LineWidth = 2;
MarkerSize = 6;

% Determine x-axis
switch data_type
    case 'N'
        % N_all = X_all;
        xlab = '$N$';
    case {'Re','p'}
        % N_all = repmat(N,1,numel(X_all));
        if strcmp(data_type,'Re')
            xlab = '$\Sigma_e$';
            X_all = Re_all;
        else
            xlab = '$p$';
            X_all = p_all;
        end        
end

% Set colours and make figure
Uf_ivs = {'iv1','iv2c','iv3c','iv4a','iv5a','iv6c'}; % visible 'groups': 3c, 1, 2c, {4a, 4c, 5a, 5c}, 6c
Yf_ivs = fieldnames(m1.Yf);
Colors = make_color_struct([Uf_ivs, Yf_ivs'], CrameriColor);

set_xlabel   = @() xlabel(xlab,'interpreter','latex','FontSize',FS_Label);
set_ylabel = @(UorY) ylabel(['$\big\| \widetilde{',UorY,'}_{\mathrm{f}}-\alpha\big\|_\mathrm{F}$'],...
    'Interpreter','latex','FontSize',FS_Label*1.25);%,'Rotation',0

fig = figure('Units',Units,'Position',FigPos);
TL0 = tiledlayout(2,1,"TileSpacing",'tight','Padding','tight');
fig.Units = 'pixels';
ax(1) = nexttile(TL0);

% --------------------- plotting for Uf_ivs -------------------------------
Uf_av_fs = replace(Uf_ivs,'iv','m');
num_Uf_ivs = numel(Uf_ivs);

% sqrt_Nall = diag(sqrt(N_all));
u_iv = struct;
u_plot_handles = [];
u_labels = {};
for kIV = 1:num_Uf_ivs
    av_name = Uf_av_fs{kIV};
    iv_name = Uf_ivs{kIV};
    
    u_iv.(av_name) = m1.Uf.(iv_name).median;

    % determine color
    switch iv_name
        case 'iv6c', label = '$\hat{\widetilde{U}}_{\mathrm{f,6c}}$'; LineStyle = '-.';  Marker = 'none';
        case 'iv2c', label = '$\hat{\widetilde{U}}_{\mathrm{f,2c}}$'; LineStyle = '-.'; Marker = 'none';
        case 'iv3c', label = '$\hat{\widetilde{U}}_{\mathrm{f,3c}}$'; LineStyle = '-';  Marker = 'none';
        case 'iv4a', label = '$\hat{\widetilde{U}}_{\mathrm{f,4a}}$'; LineStyle = '-';  Marker = 'none';
        case 'iv4c', label = '$\hat{\widetilde{U}}_{\mathrm{f,4c}}$'; LineStyle = ':';  Marker = 'o';
        case 'iv5a', label = '$\hat{\widetilde{U}}_{\mathrm{f,5a}}$'; LineStyle = ':';  Marker = 'none';
        case 'iv5c', label = '$\hat{\widetilde{U}}_{\mathrm{f,5c}}$'; LineStyle = ':';  Marker = '^';
        case 'iv1',  label = '$U_{\mathrm{f}}\vphantom{\hat{\widetilde{U}}}$'; LineStyle = '-';  Marker = 'none';
        otherwise,   label = iv_name;                                 LineStyle = '-';  Marker = 'none';
    end
    col = Colors.(iv_name);

    h = loglog(X_all, u_iv.(av_name), 'Marker', Marker, 'MarkerSize', MarkerSize,'Color', col, 'LineStyle', LineStyle, 'LineWidth', LineWidth);
    u_plot_handles = [u_plot_handles, h];
    u_labels{kIV} = label;
    hold on;
end
set_ylabel('U');
grid on;

% --------------------- plotting for Yf_ivs -------------------------------
ax(2) = nexttile(TL0);

Yf_ivs = setdiff(Yf_ivs,{'iv2b'}); % exclude case that is logically zero
Yf_av_fs = replace(Yf_ivs,'iv','m');
num_Yf_ivs = numel(Yf_ivs);

y_iv = struct;
y_plot_handles = [];
y_labels = {};
for kIV = 1:num_Yf_ivs
    av_name = Yf_av_fs{kIV};
    iv_name = Yf_ivs{kIV};
    
    y_iv.(av_name) = m1.Yf.(iv_name).median;

    % determine color
    switch iv_name
        case 'iv4b', label = '$\hat{\widetilde{Y}}_{\mathrm{f,4b}}$'; LineStyle = '-'; Marker = 'none';
        case 'iv5b', label = '$\hat{\widetilde{Y}}_{\mathrm{f,5b}}$'; LineStyle = ':'; Marker = 'none';
        case 'iv3a', label = '$\Xi_f\vphantom{\hat{\widetilde{Y}}}$'; LineStyle = '-'; Marker = 'none';
        case 'iv6a', label = '$W_f$';                                 LineStyle = '-.'; Marker = 'none';
        otherwise,   label = iv_name;                                 LineStyle = '-'; Marker = 'none';
    end
    col = Colors.(iv_name);

    h = loglog(X_all, y_iv.(av_name), 'Marker', Marker, 'MarkerSize', MarkerSize, 'Color', col, 'LineStyle', LineStyle, 'LineWidth', LineWidth);
    y_plot_handles = [y_plot_handles, h];
    y_labels{kIV} = label;
    hold on;
end
set_ylabel('Y');
set_xlabel();
grid on;
linkaxes(ax,'x');
xlim([Re_all(1) Re_all(end)]);

% fig.Units = 'points';

for k=1:2
    ax(k).FontSize = FS_Tick;
    ax(k).YLabel.FontSize = FS_Label;
end
ax(k).XLabel.FontSize = FS_Label;

%% create legends

% create 'staircase' legend for Uf_ivs
leg1 = create_columnwise_legend(ax(1), u_plot_handles, u_labels, legPosOnAx{1},...
    legendEntriesPerColumn{1} , 'Interpreter', 'latex', 'FontSize', Fs_Legend,'IconColumnWidth',15);

leg2 = make_leg2(ax(2), y_plot_handles, y_labels, legPosOnAx{2}, Fs_Legend, 2);
vertices = [leg2.Position(1:2); leg2.Position(1) leg2.Position(2)+leg2.Position(4)];
h_title = make_alpha_title_left(Fs_Legend+1,ax(2),vertices,'figure');

end

%% Local functions
function leg2 = make_leg2(ax, y_plot_handles, y_labels, leg2PosOnAx2, Fs_Legend, NumCol)
    leg2 = legend(ax, y_plot_handles, y_labels, 'Interpreter', 'latex', 'FontSize', Fs_Legend,...
    'NumColumns', NumCol,'Box', 'on','IconColumnWidth',15);
    leg2PosOnFig = RelPosAx2Fig(ax, leg2PosOnAx2, {'x','y'});
    leg2.Position([1, 2]) = leg2PosOnFig;
end

function color_struct = make_color_struct(IV_names,CrameriColor)
    nIVs = numel(IV_names);
    color_types = nan(1,nIVs);
    for k = 1:nIVs
        IV_name = IV_names{k};
        switch IV_name(3:end)
            case '1'
                col = 1;
            case {'2a','2b'}
                col = 2;
            case '2c'
                col = 3;
            case '3a'
                col = 4;
            case '3c'
                col = 5;
            case {'4a','4b'}
                col = 6;
            case '4c'
                col = 7;
            case {'5a','5b'}
                col = 8;
            case '5c'
                col = 9;
            case {'6a','6c'}
                col = 10;
            otherwise
                continue
        end
        color_types(k) = col;
    end
    color_types2 = unique(color_types(~isnan(color_types)));
    nColTypes = numel(color_types2);

    color_struct = struct;
    cCram = [crameri(CrameriColor, nColTypes-1); zeros(1,3)]; % colors
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

% Staircase Legend
function leg_handles = create_columnwise_legend(ax, handles, labels, xyPos, col_entry_counts, varargin)
    % Create a staircase-shaped legend with adjacent columns
    %
    % Input:
    %   ax                  - axes object where the legend is placed
    %   handles             - array of graphics objects to include in legend
    %   labels              - cell array of legend labels
    %   xyPos               - [x, y] position relative to axis [0,1]
    %   col_entry_counts    - vector specifying number of entries per column
    %                         e.g., [1, 2, 3] creates 3 columns with 1, 2, and 3 entries
    %   varargin            - additional legend property name-value pairs
    %
    % Legends are placed side by side, bottom-aligned. A white background patch
    % with black perimeter is drawn to create a unified staircase boundary.
    % Each legend is placed on its own invisible overlay axis.
    %
    % Usage: leg_handles = create_staircase_legend(ax, handles, labels, xyPos, col_entry_counts, ...)
    n_entries = numel(handles);
    
    % Validate col_entry_counts
    if sum(col_entry_counts) ~= n_entries
        error('Sum of col_entry_counts (%d) must equal number of legend entries (%d)', ...
              sum(col_entry_counts), n_entries);
    end
    
    % Create individual legends (one per column) with no box
    % Use overlay axes for each legend (except first which uses original axis)
    leg_handles = gobjects(numel(col_entry_counts), 1);
    leg_positions = nan(numel(col_entry_counts), 4);
    ax_overlay = gobjects(numel(col_entry_counts), 1);
    ax_overlay(1) = ax; % First legend uses original axis
    
    % Get figure and axis position (convert from layout units to figure-absolute units)
    fig = gcf;
    % Get pixel position of the axis and convert to normalized figure units
    ax_pix_pos = getpixelposition(ax);
    fig_size = getpixelposition(fig);
    % Convert pixel position to normalized figure coordinates [x, y, width, height]
    ax_pos = [ax_pix_pos(1)/fig_size(1), ax_pix_pos(2)/fig_size(2), ...
              ax_pix_pos(3)/fig_size(3), ax_pix_pos(4)/fig_size(4)];
    
    entry_idx = 1;
    
    for col = 1:numel(col_entry_counts)
        n_col_entries = col_entry_counts(col);
        col_handles = handles(entry_idx:entry_idx + n_col_entries - 1);
        col_labels = labels(entry_idx:entry_idx + n_col_entries - 1);
        
        % Create overlay axis for legends 2 onward
        if col > 1
            % Create invisible overlay axis as child of figure with same position as original
            ax_overlay(col) = axes('Parent', fig, 'Position', ax_pos,...
             'Color','none','XColor','none','YColor','none','ZColor','none',...
             'XTick',[],'YTick',[],'ZTick',[],'Box','off','HitTest','off');
            % Match limits of original axis
            ax_overlay(col).XLim = ax.XLim;
            ax_overlay(col).YLim = ax.YLim;
        end
        
        % Create single-column legend without box on the appropriate axis
        leg_handles(col) = legend(ax_overlay(col), col_handles, col_labels, 'Interpreter', 'latex',...
            'NumColumns', 1, 'Box', 'on', 'Color','white','EdgeColor','white', varargin{:});
        if col == 1
            leg_handles(col).Location = 'none';
            leg_handles(col).Position([1, 2]) = RelPosAx2Fig(ax_overlay(col), xyPos, {'x','y'});
            leg_positions(col, :) = leg_handles(col).Position;
            y_bottom = leg_positions(1, 2);
            current_x = leg_positions(1, 1) + leg_positions(1, 3); % Start to the right of first legend
        else
            old_pos = leg_handles(col).Position;
            leg_handles(col).Position = [current_x, y_bottom, old_pos(3), old_pos(4)];
            current_x = current_x + old_pos(3);
            leg_positions(col, :) = leg_handles(col).Position;
        end
        
        entry_idx = entry_idx + n_col_entries;
    end
    
    % Build staircase perimeter vertices
    function vertices = build_staircase_vertices(leg_positions)
        y_bottom = min(leg_positions(:,2));
        vertices = [];
        
        % 1. Bottom-left corner of first legend
        vertices = [vertices; leg_positions(1, 1), y_bottom];
        
        % 2. Bottom-right corner of last legend
        x_last = leg_positions(end, 1);
        w_last = leg_positions(end, 3);
        vertices = [vertices; x_last + w_last, y_bottom];
        
        % 3. Top-right corner of last legend
        h_last = leg_positions(end, 4);
        vertices = [vertices; x_last + w_last, y_bottom + h_last];
        
        % Staircase down: go left and step down for each column
        for col = size(leg_positions,1):-1:2
            x_col = leg_positions(col, 1);
            h_col = leg_positions(col, 4);
            h_prev = leg_positions(col-1, 4);
            
            % Step left to right edge of previous column
            vertices = [vertices; x_col, y_bottom + h_col];
            % Step down to top of previous column
            vertices = [vertices; x_col, y_bottom + h_prev];
        end
        
        % 4. Top-left corner of first legend
        x_first = leg_positions(1, 1);
        h_first = leg_positions(1, 4);
        vertices = [vertices; x_first, y_bottom + h_first];
    end
    vertices = build_staircase_vertices(leg_positions);
    
    % Build annotation lines for staircase perimeter
    function annoLines = build_annoLines(vertices)
        % Close the path by returning to the first vertex
        vertices = [vertices; vertices(1, :)];
        % Create staircase boundary using annotation line (with respect to axes position)
        annoLines = gobjects(size(vertices, 1)-1, 1);
        for k = 1:size(vertices,1)-1
            annoLines(k) = annotation('line', vertices(k:k+1,1), vertices(k:k+1,2), ...
                'Color', 'black', 'LineWidth', 0.5);
        end
    end
    annoLines = build_annoLines(vertices);

    % Add title annotation to the left of the staircase legend if title_left is true
    title_fontsize = leg_handles(1).FontSize + 1;
    h_title = make_alpha_title_left(title_fontsize,ax_overlay(1),[vertices(1,:);vertices(end,:)],'figure');
end