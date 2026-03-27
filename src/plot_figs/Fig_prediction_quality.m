function [fig,mean_std_ratio] = Fig_prediction_quality(fig_dir,iX,data_type,Cases,noPlotCases,insetYLims,insetPos,MainYLims,ConnectorLocation)
%   insetYLims        : 2-element cell {yLim_top, yLim_bottom} for zoom insets
%   insetPos          : 2-element cell {[X,Y,W,H]_top, [X,Y,W,H]_bottom} inset position
%                       X,Y = bottom-left corner; W,H = size; all as fraction of parent axes
%   MainYLims         : 2-element cell {yLim_top, yLim_bottom} for main axes
%                       each entry is either 'free' or a [ymin ymax] vector
%   ConnectorLocation : 2-element cell {'loc_top','loc_bottom'}, each one of:
%                       'south' - top    of rect  -> bottom of inset (inset above)
%                       'north' - bottom of rect  -> top    of inset (inset below)
%                       'west'  - right  of rect  -> left   of inset (inset right)
%                       'east'  - left   of rect  -> right  of inset (inset left)

mean_std_ratio = 0; % max. value of abs(mean)/(std. dev.)

if ~isempty(noPlotCases); Cases = setdiff(Cases,noPlotCases); end

fig_file = fullfile(fig_dir,'processed_data.mat');
load(fig_file);
load(fullfile(fig_dir,sprintf('d%s_settings.mat',data_type)));
f = opts.f;
iX_str = sprintf('iX%d',iX);

% Get number of cases
nCases = numel(Cases);

% Define base line styles and markers for variety
baseLineStyles = {'-', '--', ':', '-.'};
baseMarkers = {'none', 'o', 's', '^', 'd', 'v', '>', '<'};%, 'p', 'h', '*', 'x'};
markerSize = 5;
lineWidth = 1.2;
FS_Tick   = 15;
FS_Label  = 20;
Fs_Legend = 15;

% Get colormap from crameri - use enough colors to cycle through
nColors = max(nCases, 8); % Ensure at least 8 colors for good distribution
allColors = crameri('batlow', nColors);  % Alternative: 'roma', 'berlin', 'vik', 'oslo'

% Create cycling arrays for the actual number of cases
colors = allColors(mod(0:nCases-1, nColors) + 1, :);
lineStyles = baseLineStyles(mod(0:nCases-1, numel(baseLineStyles)) + 1);
markers = baseMarkers(mod(0:nCases-1, numel(baseMarkers)) + 1);

% Create figure with appropriate size for two-column layout
% Two-column width is typically ~7 inches
fig = figure('Units', 'inches', 'Position', [1, 1, 7, 7]);
tiledlayout(2,1,'TileSpacing','tight','Padding','tight');

% Plot mean
ax4(1) = nexttile;
hold on; box on; grid on;
h = gobjects(nCases, 1);
for kC = 1:nCases
    CaseName = Cases{kC};

    if startsWith(CaseName,'iv')
        DispName = ['IV',CaseName(3:end)];
        if strcmp(CaseName,'iv4a')
            DispName = [DispName, ' (IV2, IV4, IV5)'];
        elseif strcmp(CaseName,'iv6a')
            DispName = [DispName, ' (IV6)'];
        end
    elseif strcmp(CaseName,'CLSPC')
        DispName = 'CL-SPC';
    elseif strcmp(CaseName,'TrPred')
        DispName = 'TP';
    elseif strcmp(CaseName,'actLf')
        DispName = 'actual $L_f$';
    end
    yfhat0 = m4.yfhat.(CaseName).(iX_str).std;
    yfhat  = m4.yfhat.(CaseName).(iX_str).mean;
    h(kC) = plot(1:f, yfhat, ...
        'Color', colors(kC,:), ...
        'LineStyle', lineStyles{kC}, ...
        'Marker', markers{kC}, ...
        'MarkerSize', markerSize, ...
        'LineWidth', lineWidth, ...
        'DisplayName', DispName, ...
        'MarkerFaceColor', 'none', ...
        'MarkerEdgeColor', colors(kC,:));
    mean_std_ratio = max(mean_std_ratio,max(abs(yfhat./yfhat0),[],'all'));
end
ylabel('Mean of $\Delta\hat{y}^*_k$', 'Interpreter', 'latex', 'FontSize', FS_Label);
set(gca, 'FontSize', FS_Tick);
xlim([1, f]);

% Plot std
ax4(2) = nexttile;
hold on; box on; grid on;
for kC = 1:nCases
    CaseName = Cases{kC};
    
    if startsWith(CaseName,'iv')
        DispName = ['$j=',CaseName(3:end),'$'];
    elseif strcmp(CaseName,'CLSPC')
        DispName = 'CL-SPC';
    elseif strcmp(CaseName,'TrPred')
        DispName = 'TP';
    elseif strcmp(CaseName,'actLf')
        DispName = 'actual $L_f$';
    end
    
    yfhat2  = m4.yfhat.(CaseName).(iX_str).std;
    plot(1:f, yfhat2, ...
        'Color', colors(kC,:), ...
        'LineStyle', lineStyles{kC}, ...
        'Marker', markers{kC}, ...
        'MarkerSize', markerSize, ...
        'LineWidth', lineWidth, ...
        'DisplayName', DispName, ...
        'MarkerFaceColor', 'none', ...
        'MarkerEdgeColor', colors(kC,:));
end
xlabel('Number of time steps ahead ($k$)', 'Interpreter', 'latex', 'FontSize', FS_Label);
ylabel('Std. dev. of $\Delta\hat{y}^*_k$', 'Interpreter', 'latex', 'FontSize', FS_Label);
set(gca, 'FontSize', FS_Tick);
xlim([1, f]);

% Create a shared legend outside the plots
% Adjust number of columns based on number of cases
nLegCols = min(5, ceil(nCases/2)); % Max 5 columns, but adaptive
leg = legend(h, 'Orientation', 'horizontal', ...
    'NumColumns', nLegCols, ...
    'Location', 'northoutside', ...
    'FontSize', Fs_Legend, ...
    'Interpreter', 'latex');

% Link x-axes
linkaxes(ax4, 'x');

%% Add zoom insets

drawnow; % ensure axes positions and limits are finalised
xInset = xlim(ax4(1));
axInsets = gobjects(2,1);

for k = 1:2 % do all operations that can adjust axis positioning
    if strcmp(MainYLims{k},'free')
        MainYLims{k} = ylim(ax4(k));
    end
    if strcmp(insetYLims{k},'free')
        insetYLims{k} = ylim(ax4(k));
    end
    ax4(k).YLim = MainYLims{k}; % can adjust positioning
end
for k = 1:2
    % create flag: which axis limits are tighter
    yLimsFlag = sign(MainYLims{k}-insetYLims{k});

    % determine where to draw rectangle and connectors
    if all(yLimsFlag == [1 -1]) || all(yLimsFlag == [1 0]) || all(yLimsFlag == [0 -1])
        % main axis limits are tighter than inset limits
        % -> draw rectangle on inset, skip connectors
        drawRectangle = 'inset';
        drawConnectors = false;
    elseif all(yLimsFlag == [-1 1]) || all(yLimsFlag == [0 1]) || all(yLimsFlag == [-1 0])
        % inset axis limits are tighter than main limits
        % -> draw rectangle on main, draw connectors
        drawRectangle = 'main';
        drawConnectors = true;
    else
        % no clear containment relationship
        % -> skip rectangle and connectors entirely
        drawRectangle = 'none';
        drawConnectors = false;       
    end

    % Inset for axes
    axInsets(k) = add_zoom_inset(ax4(k), fig, FS_Tick-2, ...
        insetPos{k}(3), insetPos{k}(4), ...  % insetW, insetH  (fraction of parent axes)
        insetPos{k}(1), insetPos{k}(2), ...  % insetX, insetY  (fraction of parent axes, bottom-left corner)
        xInset, insetYLims{k});  % xInset, yInset, mainYLim

    % Drawing rectangle
    switch drawRectangle
        case 'main'
            % Draw rectangle on main
            zoomBox = rectangle(ax4(k), ...
                'Position', [xInset(1), insetYLims{k}(1), diff(xInset), diff(insetYLims{k})], ...
                'EdgeColor', [0.2 0.2 0.2], ...
                'LineStyle', '--', ...
                'LineWidth', 1.0);
        otherwise
            % skip rectangle entirely
    end

    % Draw connectors if needed
    if drawConnectors
        draw_connectors(ax4(k), axInsets(k), fig, zoomBox, ConnectorLocation{k});
    end
end

end

%% Local functions

function axInset = add_zoom_inset(axTop, fig1, insetTickFontSize, insetW, insetH, insetX, insetY, xInset, yInset)
% Define inset size/position relative to plotted area (normalized figure units)
insetW = insetW * axTop.Position(3);
insetH = insetH * axTop.Position(4);
insetX = axTop.Position(1) + insetX * axTop.Position(3);
insetY = axTop.Position(2) + insetY * axTop.Position(4);

% Create inset axes in the same figure
axInset = axes('Parent', fig1, ...
    'Units', 'normalized', ...
    'Position', [insetX insetY insetW insetH], ...
    'Box', 'on', ...
    'FontSize', insetTickFontSize);
grid on;

% Copy all plotted children from top axes into inset axes
CopiedChildren = copyobj(allchild(axTop), axInset);
axInset.ClippingStyle = 'rectangle';

% Keep copied objects that have XData with at least one non-NaN value
hasNonNanXData = false(size(CopiedChildren));
for k = 1:numel(CopiedChildren)
    if isprop(CopiedChildren(k), 'XData')
        xData = CopiedChildren(k).XData;
        if ~isempty(xData) && any(~isnan(xData(:)))
            hasNonNanXData(k) = true;
        end
    end
end
lineObjs = CopiedChildren(hasNonNanXData);
isLineOrStair = arrayfun(@(h) strcmp(get(h,'Type'),'line') || strcmp(get(h,'Type'),'stair'), lineObjs);
lineObjs = lineObjs(isLineOrStair);

for k = 1:numel(lineObjs)
    xData = lineObjs(k).XData;
    yData = lineObjs(k).YData;
    idx = find(xData >= xInset(1) & xData <= xInset(2));
    lineObjs(k).XData = xData(idx);
    lineObjs(k).YData = yData(idx);
end

% Apply final inset limits
xlim(axInset, xInset);
ylim(axInset, yInset);
end

function draw_connectors(axTop, axInset, fig1, zoomBox, connLoc)
% Draws two connector lines from a side of the zoom rectangle on the main
% axis to the corresponding side of the inset axis.
%
% connLoc : one of 'south','north','west','east'
%   'south' - top    of rectangle -> bottom of inset  (inset above)
%   'north' - bottom of rectangle -> top    of inset  (inset below)
%   'west'  - right  of rectangle -> left   of inset  (inset to the right)
%   'east'  - left   of rectangle -> right  of inset  (inset to the left)

if nargin < 5 || isempty(connLoc)
    connLoc = 'south';
end

% Move ticks to opposite side of connection point
switch lower(connLoc)
    case 'south', axInset.XAxisLocation = 'top';
    case 'north', axInset.XAxisLocation = 'bottom';
    case 'west',  axInset.YAxisLocation = 'right';
    case 'east',  axInset.YAxisLocation = 'left';
end

% Create two annotation lines (initial dummy positions; updated immediately)
hLine1 = annotation(fig1, 'line', [0 0], [0 0], ...
    'Color', zoomBox.EdgeColor, 'LineWidth', 0.9);
hLine2 = annotation(fig1, 'line', [0 0], [0 0], ...
    'Color', zoomBox.EdgeColor, 'LineWidth', 0.9);

    function update_lines(~,~)
        [p1a, p1b, p2a, p2b] = compute_connector_pts(axTop, axInset, zoomBox, connLoc);
        hLine1.X = [p1a(1), p1b(1)];
        hLine1.Y = [p1a(2), p1b(2)];
        hLine2.X = [p2a(1), p2b(1)];
        hLine2.Y = [p2a(2), p2b(2)];
    end

% Initial draw
update_lines();

% Listen to position changes on axes and figure resize so lines stay in sync
addlistener(axTop,   'Position', 'PostSet', @update_lines);
addlistener(axInset, 'Position', 'PostSet', @update_lines);
addlistener(fig1, 'SizeChanged', @update_lines);
end

function [p1a, p1b, p2a, p2b] = compute_connector_pts(axTop, axInset, zoomBox, connLoc)
% Returns two pairs of points (in normalised figure coords) that form
% the two connector lines between the zoom rectangle and the inset axes.

% Derive zoom region from rectangle position
xInset = [zoomBox.Position(1), zoomBox.Position(1) + zoomBox.Position(3)];
yInset = [zoomBox.Position(2), zoomBox.Position(2) + zoomBox.Position(4)];

xTopL    = xlim(axTop);
yTopL    = ylim(axTop);
axTopPos = axTop.Position;
axInsPos = axInset.Position;

% Rectangle corners in normalised figure coordinates
x1n = axTopPos(1) + (xInset(1)-xTopL(1))/diff(xTopL) * axTopPos(3);
x2n = axTopPos(1) + (xInset(2)-xTopL(1))/diff(xTopL) * axTopPos(3);
if strcmpi(axTop.YDir, 'reverse')
    y1n = axTopPos(2) + (yTopL(2)-yInset(1))/diff(yTopL) * axTopPos(4);
    y2n = axTopPos(2) + (yTopL(2)-yInset(2))/diff(yTopL) * axTopPos(4);
else
    y1n = axTopPos(2) + (yInset(1)-yTopL(1))/diff(yTopL) * axTopPos(4);
    y2n = axTopPos(2) + (yInset(2)-yTopL(1))/diff(yTopL) * axTopPos(4);
end
yBot = min(y1n, y2n);
yTop = max(y1n, y2n);

% Rectangle side endpoints
rectTopLeft     = [x1n, yTop];
rectTopRight    = [x2n, yTop];
rectBottomLeft  = [x1n, yBot];
rectBottomRight = [x2n, yBot];
rectRightTop    = [x2n, yTop];
rectRightBottom = [x2n, yBot];
rectLeftTop     = [x1n, yTop];
rectLeftBottom  = [x1n, yBot];

% Inset side endpoints in normalised figure coordinates
iL = axInsPos(1);               iR = axInsPos(1) + axInsPos(3);
iB = axInsPos(2);               iT = axInsPos(2) + axInsPos(4);
insetBotLeft   = [iL, iB];      insetBotRight  = [iR, iB];
insetTopLeft   = [iL, iT];      insetTopRight  = [iR, iT];

switch lower(connLoc)
    case 'south'   % top of rect -> bottom of inset
        boxPts   = [rectTopLeft;    rectTopRight];
        insetPts = [insetBotLeft;   insetBotRight];
    case 'north'   % bottom of rect -> top of inset
        boxPts   = [rectBottomLeft; rectBottomRight];
        insetPts = [insetTopLeft;   insetTopRight];
    case 'west'    % right of rect -> left of inset
        boxPts   = [rectRightTop;   rectRightBottom];
        insetPts = [insetTopLeft;   insetBotLeft];
    case 'east'    % left of rect -> right of inset
        boxPts   = [rectLeftTop;    rectLeftBottom];
        insetPts = [insetTopRight;  insetBotRight];
    otherwise
        error('ConnectorLocation must be ''south'', ''north'', ''west'', or ''east''.');
end

% Pair endpoints to avoid crossing connectors
d11 = norm(boxPts(1,:) - insetPts(1,:));
d22 = norm(boxPts(2,:) - insetPts(2,:));
d12 = norm(boxPts(1,:) - insetPts(2,:));
d21 = norm(boxPts(2,:) - insetPts(1,:));
if (d11 + d22) <= (d12 + d21)
    p1a = boxPts(1,:);  p1b = insetPts(1,:);
    p2a = boxPts(2,:);  p2b = insetPts(2,:);
else
    p1a = boxPts(1,:);  p1b = insetPts(2,:);
    p2a = boxPts(2,:);  p2b = insetPts(1,:);
end
end