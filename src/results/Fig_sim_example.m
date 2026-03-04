function fig1 = Fig_sim_example(data_type,iX,ks,noPlotCases,subdir1,marker_interval)

plot_run; % create plot of Monte-Carlo simulation

% 
lgd = findobj(gcf, 'Type', 'Legend');
lgd.NumColumns = 3;
lgd.Location = 'southwest';
lgd.FontSize = 12;

fig1 = gcf;
fig1.Position = [680   283   940   830];%595];

allAxes = findall(gcf, 'Type', 'axes');

axisFontSize = 14;
% Loop through and set label font sizes
for i = 1:length(allAxes)
    allAxes(i).FontSize = axisFontSize-2; % numbers on axis
    allAxes(i).XLabel.FontSize = axisFontSize;
    allAxes(i).YLabel.FontSize = axisFontSize;
end

xlim(gca,[0 length(yr0)+length(yr1)-1])

%% inset axes
% Select the top-most subplot (largest y-position) as source for zoom-in
plotAxes = allAxes(isgraphics(allAxes, 'axes'));
axPos = reshape([plotAxes.Position], 4, []).';
[~, idxTop] = max(axPos(:,2)); % find top axes
[~, idxBot] = min(axPos(:,2)); % find bottom axes
axTop = plotAxes(idxTop);
axBot = plotAxes(idxBot);
axBot.YLim = [-16 16];

% Define inset size/position relative to axes (normalized figure units)
insetW = 0.32; % this is the width of the inset as a fraction of the axes width
insetH = 0.3;  % this is the height of the inset as a fraction of the axes height
insetX = 0.05; % distance from left edge of axes; this is added to axes left pos. below
insetY = 0.68;  % distance from bottom edge of axes; this is added to axes bottom pos. below

xInset = [1210 1370];
yInsetTop = [9 12.6];
yInsetBot = [8.5 12.2];

add_zoom_inset(axTop, fig1, axisFontSize, insetW, insetH, insetX, insetY, xInset, yInsetTop);
add_zoom_inset(axBot, fig1, axisFontSize, insetW, insetH, insetX, insetY, xInset, yInsetBot);

end

function add_zoom_inset(axTop, fig1, axisFontSize, insetW, insetH, insetX, insetY, xInset, yInset)
% Define inset size/position relative to axes (normalized figure units)
insetW = insetW * axTop.Position(3);
insetH = insetH * axTop.Position(4);
insetX = axTop.Position(1) + insetX * axTop.Position(3);
insetY = axTop.Position(2) + insetY * axTop.Position(4);

% Create inset axes in the same figure
axInset = axes('Parent', fig1, ...
    'Position', [insetX insetY insetW insetH], ...
    'Box', 'on', ...
    'FontSize', axisFontSize-3);

% Copy all plotted children from top axes into inset axes
CopiedChildren = copyobj(allchild(axTop), axInset);
axInset.XGrid = axTop.XGrid;
axInset.YGrid = axTop.YGrid;
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

% Zoom inset to the last 20% of the x-range of the top axes
xTop = xlim(axTop);
yTop = ylim(axTop);
if nargin < 8 || isempty(xInset)
    xSpan = diff(xTop);
    xInset = [xTop(2) - 0.20*xSpan, xTop(2)];
    if xInset(1) >= xInset(2)
        xInset = xTop;
    end
end

for k = 1:numel(lineObjs)
    xData = lineObjs(k).XData;
    yData = lineObjs(k).YData;
    idx = find(xData >= xInset(1) & xData <= xInset(2));
    lineObjs(k).XData = xData(idx);
    lineObjs(k).YData = yData(idx);
end

if nargin < 9 || isempty(yInset)
    % Compute y-limits in the inset based only on data within xInset
    yMin = inf;
    yMax = -inf;
    for k = 1:numel(lineObjs)
        yData = lineObjs(k).YData;
        yData = yData(~isnan(yData)); % ignore NaNs
        if ~isempty(yData)
            yMin = min(yMin, min(yData));
            yMax = max(yMax, max(yData));
        end
    end

    % Add a small padding around y-data; otherwise fall back to top-axes y-limits
    if isfinite(yMin) && isfinite(yMax) && yMax > yMin
        yPad = 0.08*(yMax - yMin);
        yInset = [yMin - yPad, yMax + yPad];
    else
        yInset = yTop;
    end
end

% Apply final inset limits
xlim(axInset, xInset);
ylim(axInset, yInset);

% Draw a box on the top axes indicating the zoomed region
zoomBox = rectangle(axTop, ...
    'Position', [xInset(1), yInset(1), diff(xInset), diff(yInset)], ...
    'EdgeColor', [0.2 0.2 0.2], ...
    'LineStyle', '--', ...
    'LineWidth', 1.0);

% Compute box corners in normalized figure coordinates
xTopL = xlim(axTop);
yTopL = ylim(axTop);
axTopPos = axTop.Position;
axInsetPos = axInset.Position;

x1n = axTopPos(1) + (xInset(1)-xTopL(1))/diff(xTopL) * axTopPos(3);
x2n = axTopPos(1) + (xInset(2)-xTopL(1))/diff(xTopL) * axTopPos(3);
if strcmpi(axTop.YDir, 'reverse')
    y1n = axTopPos(2) + (yTopL(2)-yInset(1))/diff(yTopL) * axTopPos(4);
    y2n = axTopPos(2) + (yTopL(2)-yInset(2))/diff(yTopL) * axTopPos(4);
else
    y1n = axTopPos(2) + (yInset(1)-yTopL(1))/diff(yTopL) * axTopPos(4);
    y2n = axTopPos(2) + (yInset(2)-yTopL(1))/diff(yTopL) * axTopPos(4);
end

boxLeftBottom  = [x1n min(y1n,y2n)];
boxLeftTop     = [x1n max(y1n,y2n)];
boxRightBottom = [x2n min(y1n,y2n)];
boxRightTop    = [x2n max(y1n,y2n)];
boxBottomLeft  = boxLeftBottom;
boxBottomRight = boxRightBottom;
boxTopLeft     = boxLeftTop;
boxTopRight    = boxRightTop;

% Inset corners in normalized figure coordinates
insetLeftBottom  = [axInsetPos(1), axInsetPos(2)];
insetLeftTop     = [axInsetPos(1), axInsetPos(2)+axInsetPos(4)];
insetRightBottom = [axInsetPos(1)+axInsetPos(3), axInsetPos(2)];
insetRightTop    = [axInsetPos(1)+axInsetPos(3), axInsetPos(2)+axInsetPos(4)];

% Select which sides to connect based on relative inset position
boxCenter = [mean([x1n x2n]), mean([min(y1n,y2n) max(y1n,y2n)])];
insetCenter = [axInsetPos(1)+0.5*axInsetPos(3), axInsetPos(2)+0.5*axInsetPos(4)];
delta = insetCenter - boxCenter;

if abs(delta(1)) >= abs(delta(2))
    if delta(1) >= 0
        boxPts = [boxRightTop; boxRightBottom];
        insetPts = [insetLeftTop; insetLeftBottom];
    else
        boxPts = [boxLeftTop; boxLeftBottom];
        insetPts = [insetRightTop; insetRightBottom];
    end
else
    if delta(2) >= 0
        boxPts = [boxTopLeft; boxTopRight];
        insetPts = [insetLeftBottom; insetRightBottom];
    else
        boxPts = [boxBottomLeft; boxBottomRight];
        insetPts = [insetLeftTop; insetRightTop];
    end
end

% Pair endpoints to avoid crossing connectors
d11 = norm(boxPts(1,:) - insetPts(1,:));
d22 = norm(boxPts(2,:) - insetPts(2,:));
d12 = norm(boxPts(1,:) - insetPts(2,:));
d21 = norm(boxPts(2,:) - insetPts(1,:));
if (d11 + d22) <= (d12 + d21)
    p1a = boxPts(1,:); p1b = insetPts(1,:);
    p2a = boxPts(2,:); p2b = insetPts(2,:);
else
    p1a = boxPts(1,:); p1b = insetPts(2,:);
    p2a = boxPts(2,:); p2b = insetPts(1,:);
end

annotation(fig1, 'line', [p1a(1) p1b(1)], [p1a(2) p1b(2)], ...
    'Color', zoomBox.EdgeColor, 'LineWidth', 0.9);
annotation(fig1, 'line', [p2a(1) p2b(1)], [p2a(2) p2b(2)], ...
    'Color', zoomBox.EdgeColor, 'LineWidth', 0.9);
end