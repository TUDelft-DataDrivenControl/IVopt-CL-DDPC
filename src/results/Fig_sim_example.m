function fig1 = Fig_sim_example(data_type,iX,ks,noPlotCases,subdir1,marker_interval)

FS_Tick   = 8;
FS_Label  = 9;
Fs_Legend = 7;

%% Choose settings for which to show simulation results
%% navigate to subdir1
if isempty(subdir1)
    [subdir1,~] = get_subdir1(data_type);
end
pwd1 = pwd;
cd(subdir1);

%% load data
load(sprintf('d%s_settings.mat',data_type));

subdir2s = dir(subdir1);
isub = [subdir2s(:).isdir]; 
subdir2s = {subdir2s(isub).name};
subdir2s = subdir2s(~ismember(subdir2s,{'.','..','mfiles'}));
switch data_type
    case {'N','p'}
        if strcmp(data_type,'N')
            N = N_all(iX); X = N;
        else % -> 'p'
            p = p_all(iX); X = p;
        end
        subdir2 = choose_subdir_by_number(subdir2s, X);
        cd(subdir2); % navigate into <subdir2>

        % load file with settings in <subdir2>
        settingsFile = find_settingsFile();
        load(settingsFile);
    
    case 'Re'
        subdir2 = choose_subdir_by_iX(subdir2s, iX);
        cd(subdir2); % navigate into <subdir2>
end
load(sprintf('seed_%d.mat',seeds(ks,iX)));

f = opts.f;
Ncl = opts.Ncl;

%% Plotting

% determine unstable cases, do not plot those
fns = fieldnames(Tcl);
for k = 1:numel(fns)
    fn = fns{k};
    if ~isstable(Tcl.(fn))
        noPlotCases = [noPlotCases,{fn}];
    end
end
noPlotCases = unique(noPlotCases);

% select cases to plot
Cases  = setdiff(Cases,noPlotCases);
nCases = numel(Cases);

% determine colours to use
cCram = crameri('batlow', nCases); % colors

% time steps
Tsteps0 = 0:Nbar-1;
Tsteps1 = Nbar:Nbar+Ncl-1;
Tsteps  = [Tsteps0 Tsteps1];

% make figure and tiledlayout
fig1 = figure();
fig1.OuterPosition(3:4) = [700 385];
nColsPerTile = 1;
nTileCols = nColsPerTile*3;
nTileRows = 2;
tl = tiledlayout(nTileRows,nTileCols,'TileSpacing','compact','Padding','tight');

% Define line styles and marker styles to cycle through
line_styles = {'-'};%, '--', '-.'};
marker_styles = {'o', 's', 'd', '^', 'none','v', '>', '<', 'p', 'h'};
% ------------------------------- outputs ----------------------------
ax11 = nexttile(1,[1 2*nColsPerTile]);
stairs(Tsteps,[e0 e1],'r-','DisplayName','$e_k$'); hold on; % innovation noise
stairs(Tsteps0,y0,'b-','DisplayName','past data');
plot_uycl(Cases,true,line_styles,marker_styles,cCram,Tsteps1,y_cl,marker_interval);
stairs(0:Nbar+Ncl-1+f,[yr0 yr1],'k-','LineWidth',1.5,'DisplayName','ref.');  % references

xline(Nbar-0.5,'LineWidth',2,'Color','k','LineStyle','-.','HandleVisibility','off');
grid on;
ylabel('$y_k$','Interpreter','latex','FontSize',12);

% add a legend north of the plots
lgd = legend('Interpreter','latex','Orientation','Horizontal',...
    'NumColumns',8,'FontSize',Fs_Legend,'IconColumnWidth',20);
lgd.Layout.Tile = 'north';

% focus box:
ax11r = nexttile(2*nColsPerTile+1,[1 nColsPerTile]);
AxCopy1 = copyobj(allchild(ax11), ax11r);
grid on; box on;
xInset = [2210 2370]; yInsetTop = [8.9 13.4];%[1210 1370]
xlim(ax11r, xInset); ylim(ax11r, yInsetTop);
ax11r.YAxisLocation = 'right';
trim_copied_objects(AxCopy1, xInset, yInsetTop);

% Draw a box on indicating the zoomed region
zoomBox = rectangle(ax11, ...
    'Position', [xInset(1), yInsetTop(1), diff(xInset), diff(yInsetTop)], ...
    'EdgeColor', [0.5 0.5 0.5], ...
    'LineStyle', '-', ...
    'LineWidth', 1.0);

% ------------------------------- inputs -----------------------------
ax12 = nexttile(nTileCols+1,[1 2*nColsPerTile]);
stairs(Tsteps0,u0,'b-'); hold on;
plot_uycl(Cases,false,line_styles,marker_styles,cCram,Tsteps1,u_cl,marker_interval);
stairs(Nbar:Nbar+Ncl-1+f,ur1,'k-','LineWidth',1.5);  % references

xline(Nbar-0.5,'LineWidth',2,'Color','k','LineStyle','-.');
grid on;
ylabel('$u_k$','Interpreter','latex','FontSize',12);
xlabel('Time step','Interpreter','latex','FontSize',12);

linkaxes([ax11 ax12],'x');
xlim(ax12,[floor(Nbar*0.8) length(yr0)+length(yr1)-1]);
ylim(ax12,[-16 14.5]);

% focus box:
ax12r = nexttile(nTileCols+2*nColsPerTile+1,[1 nColsPerTile]);
AxCopy2 = copyobj(allchild(ax12), ax12r);
grid on; box on
yInsetBot = [8.5 12.6];
xlim(ax12r, xInset); ylim(ax12r, yInsetBot);
ax12r.YAxisLocation = 'right';
trim_copied_objects(AxCopy2, xInset, yInsetBot);

% Draw a box on indicating the zoomed region
zoomBox = rectangle(ax12, ...
    'Position', [xInset(1), yInsetBot(1), diff(xInset), diff(yInsetBot)], ...
    'EdgeColor', [0.5 0.5 0.5], ...
    'LineStyle', '-', ...
    'LineWidth', 1.0);

%% rescaling figure
axRelHeight = [ax11.Position(4) ax12.Position(4)]; % store old relative axis height
% resize figure to be 1 column width (252 pt)
fig1.Units = 'points';
ptwidth = 516; % 252 for single-column
fig1.OuterPosition(3:4) = [ptwidth, ptwidth/fig1.OuterPosition(3)*fig1.OuterPosition(4)];
axHeightScaleFac = max(axRelHeight./[ax11.Position(4) ax12.Position(4)]);
fig1.OuterPosition(4) = fig1.OuterPosition(4)*axHeightScaleFac; % restore old relative height of axis after resizing figure

allAxes = findall(fig1, 'Type', 'axes');
% Loop through and set label font sizes
for i = 1:length(allAxes)
    allAxes(i).FontSize = FS_Tick; % numbers on axis
    allAxes(i).XLabel.FontSize = FS_Label;
    allAxes(i).YLabel.FontSize = FS_Label;
end

%% connector lines from zoom boxes to zoom tile axes
drawnow;  % ensure axis Positions are current before computing figure coordinates
draw_zoom_connectors(fig1, ax11, ax11r, xInset, yInsetTop);
draw_zoom_connectors(fig1, ax12, ax12r, xInset, yInsetBot);

%% relocate legend
tl.PositionConstraint = 'innerposition';
lgd.Location = 'none';
lgd.Position(2) = lgd.Position(2)-0.07;
%% inset axes
% Select the top-most subplot (largest y-position) as source for zoom-in
% plotAxes = allAxes(isgraphics(allAxes, 'axes'));
% axPos = reshape([plotAxes.Position], 4, []).';
% [~, idxTop] = max(axPos(:,2)); % find top axes
% [~, idxBot] = min(axPos(:,2)); % find bottom axes
% axTop = plotAxes(idxTop);
% axBot = plotAxes(idxBot);
% axBot.YLim = [-17 17];
% axBot.YLim = axBot.YLim+3;
% axTop.YLim = axTop.YLim+3;
% 
% % Define inset size/position relative to axes (normalized figure units)
% insetW = 0.3; % this is the width of the inset as a fraction of the axes width
% insetH = 0.8;  % this is the height of the inset as a fraction of the axes height
% insetX = 0.075; % distance from left edge of axes; this is added to axes left pos. below
% insetY = 0.15;  % distance from bottom edge of axes; this is added to axes bottom pos. below
% 
% xInset = [1210 1370];
% yInsetTop = [9 12.6];
% yInsetBot = [8.5 12.2];
% pause(0.5)
% 
% add_zoom_inset(axTop, fig1, FS_Tick, insetW, insetH, insetX, insetY, xInset, yInsetTop);
% add_zoom_inset(axBot, fig1, FS_Tick, insetW, insetH, insetX, insetY, xInset, yInsetBot);

cd(pwd1);
end
%% Helper functions
function plot_uycl(Cases,legend_flag,line_styles,marker_styles,cCram,Tsteps1,y_cl,marker_interval)
    nCases = numel(Cases);
    for k = 1:nCases
        CaseName = Cases{k};
                
        % Cycle through line & marker styles, & marker indices
        line_style = line_styles{mod(k-1, numel(line_styles)) + 1};
        marker_style = marker_styles{mod(k-1, numel(marker_styles)) + 1};
        marker_indices = get_marker_indices(k, nCases, marker_interval, numel(Tsteps1));

        % Add markers at specified intervals
        hold on;
        if legend_flag
            switch CaseName
                case 'CLSPC'
                    DispName = 'CL-SPC';
                case 'TrPred'
                    DispName = 'TP';
                case 'actLf'
                    DispName = 'actual $L_f$';
                otherwise
                    if startsWith(CaseName,'iv')
                        DispName = ['IV',CaseName(3:end)];
                    else
                        DispName = CaseName;
                    end
            end
            plot(nan, nan, ...      just for legend
                'Color', cCram(k,:), ...
                'LineStyle', line_style, ...
                'Marker', marker_style, ...
                'MarkerSize', 3,...
                'DisplayName',DispName);
        end
        stairs(Tsteps1, y_cl.(CaseName), ... plotting of stairs
            'Color', cCram(k,:), ...
            'LineStyle', line_style, ...
            'HandleVisibility', 'off');

        plot(Tsteps1(marker_indices), y_cl.(CaseName)(marker_indices), ...
            'Color', cCram(k,:), ...
            'LineStyle', 'none', ...
            'Marker', marker_style, ...
            'MarkerSize', 3,...
            'HandleVisibility', 'off');
    end
end

function add_zoom_inset(axTop, fig1, insetTickFontSize, insetW, insetH, insetX, insetY, xInset, yInset)
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
% axInset.XGrid = axTop.XGrid;
% axInset.YGrid = axTop.YGrid;
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

function draw_zoom_connectors(fig, axSrc, axZoom, xInset, yInset)
% DRAW_ZOOM_CONNECTORS  Draw two annotation lines from the right edge of
% the zoom rectangle (on axSrc, in data coordinates xInset/yInset) to the
% left edge of the zoom tile axis (axZoom). Connector colour matches the
% zoom-box rectangle already drawn on axSrc.
%
% Both axes must already have their final Position set (call drawnow first).
% Works correctly regardless of the TiledChartLayout PositionConstraint.

% Retrieve the existing zoom-box rectangle on axSrc (last rectangle added)
rects = findobj(axSrc, 'Type', 'rectangle');
if isempty(rects)
    lineColor = [0.5 0.5 0.5];
else
    lineColor = rects(1).EdgeColor;
end

xSrcLim = xlim(axSrc);
ySrcLim = ylim(axSrc);

% InnerPosition is always the plot-box in normalised figure units,
% regardless of PositionConstraint ('innerposition' or 'outerposition').
oldUnits = axSrc.Units;
axSrc.Units = 'normalized';
axSrcInner = axSrc.InnerPosition;   % [x y w h] of plot area in figure units
axSrc.Units = oldUnits;

% Right edge of zoom box -> normalised figure x (using the plot-area width)
x_right = axSrcInner(1) + (xInset(2) - xSrcLim(1)) / diff(xSrcLim) * axSrcInner(3);

% Top and bottom of zoom box -> normalised figure y
if strcmpi(axSrc.YDir, 'reverse')
    y_top = axSrcInner(2) + (ySrcLim(2) - yInset(1)) / diff(ySrcLim) * axSrcInner(4);
    y_bot = axSrcInner(2) + (ySrcLim(2) - yInset(2)) / diff(ySrcLim) * axSrcInner(4);
else
    y_top = axSrcInner(2) + (yInset(2) - ySrcLim(1)) / diff(ySrcLim) * axSrcInner(4);
    y_bot = axSrcInner(2) + (yInset(1) - ySrcLim(1)) / diff(ySrcLim) * axSrcInner(4);
end

% Left edge of zoom tile -> normalised figure coordinates (OuterPosition
% gives the tile bounding box; InnerPosition gives the plot area).
oldUnits = axZoom.Units;
axZoom.Units = 'normalized';
axZoomInner = axZoom.InnerPosition;
axZoom.Units = oldUnits;

x_left = axZoomInner(1);
z_top  = axZoomInner(2) + axZoomInner(4);
z_bot  = axZoomInner(2);

annotation(fig, 'line', [x_right x_left], [y_top z_top], ...
    'Color', lineColor, 'LineWidth', 0.9);
annotation(fig, 'line', [x_right x_left], [y_bot z_bot], ...
    'Color', lineColor, 'LineWidth', 0.9);
end

function trim_copied_objects(objs, xLim, yLim)
% TRIM_COPIED_OBJECTS  Remove graphics objects that lie entirely outside
% [xLim(1) xLim(2)] x [yLim(1) yLim(2)], and trim the XData/YData of
% remaining line/stair objects so only the data within the visible range
% (plus one boundary neighbour on each side) is kept.
for k = numel(objs):-1:1   % reverse so deletions don't shift indices
    obj = objs(k);
    objType = get(obj, 'Type');

    % --- Constant vertical lines (xline) ---
    if strcmpi(objType, 'constantline')
        if obj.Value < xLim(1) || obj.Value > xLim(2)
            delete(obj);
        end
        continue;
    end

    % --- Objects without XData / YData (patches, rectangles, text, ...) ---
    % MATLAB clips these at the axis limits automatically; keep them.
    if ~isprop(obj,'XData') || ~isprop(obj,'YData')
        continue;
    end

    xData = obj.XData(:)';
    yData = obj.YData(:)';

    % Delete empty or all-NaN objects
    if isempty(xData) || all(isnan(xData))
        delete(obj);
        continue;
    end

    % Delete objects whose bounding box does not overlap the visible region
    if max(xData) < xLim(1) || min(xData) > xLim(2) || ...
       max(yData) < yLim(1) || min(yData) > yLim(2)
        delete(obj);
        continue;
    end

    % Trim to x range; keep one neighbouring point on each side so that
    % line / stair segments that cross the boundary remain intact.
    inX    = xData >= xLim(1) & xData <= xLim(2);
    iFirst = find(inX, 1, 'first');
    iLast  = find(inX, 1, 'last');
    if isempty(iFirst)
        delete(obj);
        continue;
    end
    iStart = max(1, iFirst - 1);
    iEnd   = min(numel(xData), iLast + 1);
    obj.XData = xData(iStart:iEnd);
    obj.YData = yData(iStart:iEnd);
end
end

function marker_indices = get_marker_indices(caseIdx, nCases, marker_interval, nPoints)
startIdx = floor((caseIdx-1)/nCases*marker_interval) + 1;
marker_indices = startIdx:marker_interval:nPoints;
if isempty(marker_indices)
    marker_indices = 1;
else
    marker_indices = marker_indices(marker_indices <= nPoints);
end
end

function chosenDir = choose_subdir_by_number(subdirs, targetNum)
% CHOOSE_SUBDIR_BY_NUMBER selects the subdirectory whose trailing integer
% matches targetNum.
%
% Inputs:
%   subdirs   - cell array of subdirectory names (strings)
%   targetNum - integer to match at the end of the subdirectory name
%
% Output:
%   chosenDir - matching subdirectory name (string). Empty if no match.

    chosenDir = '';  % default (if no match)

    for i = 1:numel(subdirs)
        % Extract trailing number using regexp
        tokens = regexp(subdirs{i}, '(\d+)$', 'tokens');
        if ~isempty(tokens)
            num = str2double(tokens{1}{1});
            if num == targetNum
                chosenDir = subdirs{i};
                return;  % stop after first match
            end
        end
    end

    if isempty(chosenDir)
        warning('No subdirectory ends with the number %d.', targetNum);
    end
end

function settingsFile = find_settingsFile(dirPath)
% FIND_MAT_FILES returns all .mat files in the specified directory ending
% with 'settings.mat'.
%
% Inputs:
%   dirPath - path to the directory (string). Defaults to pwd.
%
% Outputs:
%   settingsFiles - cell array of full paths to .mat files ending with 'settings.mat'

    if nargin < 1
        dirPath = pwd;  % default to current directory
    end

    % Get all .mat files
    settingsFile = dir(fullfile(dirPath, '*settings.mat'));

    % convert to cell array
    settingsFile = {settingsFile.name};

    if numel(settingsFile) > 1
        error('only expecting one settings file')
    end
    settingsFile = settingsFile{1};
end

function chosenDir = choose_subdir_by_iX(subdirs, iX)
% CHOOSE_SUBDIR_BY_IX selects the subdirectory whose leading zero-padded number matches iX.
% Example: for iX=3, matches '003_someName' if present.
    chosenDir = '';
    for i = 1:numel(subdirs)
        dirName = subdirs{i};
        tokens = regexp(dirName, '^(\d+)', 'tokens');
        if ~isempty(tokens)
            dirNum = str2double(tokens{1}{1});
            if dirNum == iX
                chosenDir = dirName;
                return;
            end
        end
    end
    if isempty(chosenDir)
        warning('No subdirectory found starting with number %d.', iX);
    end
end