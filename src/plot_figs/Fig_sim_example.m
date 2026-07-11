function fig1 = Fig_sim_example(data_type,iX,ks,noPlotCases,subdir1,marker_interval,ylim_y1,ylim_y2)

FS_Tick   = 8;
FS_Label  = 9;
Fs_Legend = 7;

%% navigate to subdir1
pwd1 = pwd;
cd(subdir1);

%% load data
load(sprintf('d%s_settings.mat',data_type));

subdir2s = dir(subdir1);
isub = [subdir2s(:).isdir]; 
subdir2s = {subdir2s(isub).name};
subdir2s = subdir2s(~ismember(subdir2s,{'.','..'}));
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
Cases  = setdiff(opts.Cases,noPlotCases);
nCases = numel(Cases);

% determine colours to use
cCram = crameri('batlow', nCases); % colors
cCram(end,:) = cCram(end,:)/2;     % to improve visibility

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
ylim(ylim_y1{1});

% add a legend north of the plots
lgd = legend('Interpreter','latex','Orientation','Horizontal',...
    'NumColumns',8,'FontSize',Fs_Legend,'IconColumnWidth',20);
lgd.Layout.Tile = 'north';

% focus box:
ax11r = nexttile(2*nColsPerTile+1,[1 nColsPerTile]);
AxCopy1 = copyobj(allchild(ax11), ax11r);
grid on; box on;
xInset = [2210 2370] + (opts.N-1e3 + opts.p-20 + opts.f-20); yInsetTop = ylim_y1{2};
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
ylim(ax12,ylim_y2{1});

% focus box:
ax12r = nexttile(nTileCols+2*nColsPerTile+1,[1 nColsPerTile]);
AxCopy2 = copyobj(allchild(ax12), ax12r);
grid on; box on
yInsetBot = ylim_y2{2};
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

cd(pwd1);
end
%% Local functions
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

function draw_zoom_connectors(fig, axSrc, axZoom, xInset, yInset)
% DRAW_ZOOM_CONNECTORS  Draw two annotation lines from the right edge of
% the zoom rectangle (on axSrc, in data coordinates xInset/yInset) to the
% left edge of the zoom tile axis (axZoom). Connector colour matches the
% zoom-box rectangle already drawn on axSrc.
%
% Both axes must already have their final Position set (call drawnow first).

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