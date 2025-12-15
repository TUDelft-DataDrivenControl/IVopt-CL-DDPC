function axLeg = customLegend(entries, hParentAx, opts)
%CUSTOMLEGEND Create a custom legend with patch+line+LaTeX text
% The legend auto-resizes/repositions itself to fit all text inside.
%
%   axLeg = customLegend(parentAx, entries)
%   axLeg = customLegend(parentAx, entries, Name,Value,...)
%
% Inputs:
%   parentAx  : axes handle
%   entries   : struct array with fields:
%                 .Color      : RGB triplet or short color name
%                 .Alpha      : Face transparency (0–1)
%                 .LineStyle  : Line style string (e.g. '-', '--', ':')
%                 .LineWidth  : Line width in points
%                 .Text       : Label string (LaTeX interpreted)
%
% Name-Value options (all optional):
%   'Position'        : [x y w h] position in normalized figure units.
%                       If specified, this overrides 'Location'.
%   'Location'        : Legend placement keyword:
%                       "northwest","northeast","southwest","southeast",
%                       "north","south","east","west","center".
%                       Default is "northeast", matching MATLAB's legend.
%   'XPadding'        : Horizontal padding (default 0.05)
%   'YPadding'        : Vertical padding (default 0.2)
%   'BoxWidth'        : Width of color patch (default 0.3)
%   'RowHeightFactor' : Row height relative to max text height (default 1)
%   'BoxHeightFactor' : Patch height relative to max text height (default 1)
%   'YSpacingFactor'  : Vertical spacing relative to max text height (default 0.2)
%   'FontSize'        : Font size for labels (default 12)
%
% Output:
%   axLeg : handle to the legend axes

%% Argument validation
arguments
    entries (1,:) struct
    hParentAx {mustBeA(hParentAx,"matlab.graphics.axis.Axes")} = gca
    opts.Position = []; % explicit [x y w h] relative to axis, empty = auto
    opts.XPadding (1,1) double = 5  % [pixels]
    opts.YPadding (1,1) double = 5  % [pixels]
    opts.BoxWidth (1,1) double = 30 % [pixels]
    opts.RowHeightFactor (1,1) double = 1
    opts.BoxHeightFactor (1,1) double = 1
    opts.YSpacingFactor (1,1) double = 0.2
    opts.MaxRelWidth  double = []; % max relative width "of legend relative to axes" (default = 0.2 * opts.cols)
    opts.MaxRelHeight double = []; % max relative height ""                          (default = 1)
    opts.RelScaling (1,1) logical = true;        % apply above relative scaling to potentially enlarge figure
    opts.FontSize (1,1) double = 12
    opts.Location (1,1) string {mustBeMember(opts.Location, ...
        ["northwest","northeast","southwest","southeast","north","south","east","west","center"])} = "northeast"
    opts.cols (1,1) double {mustBeInteger,mustBeGreaterThanOrEqual(opts.cols,1)} = 1
end
opts = validateOptions(opts,entries);

%% Preliminaries
pos0      = opts.Position;
fontSize  = opts.FontSize;
xPad_pix  = opts.XPadding;
yPad_pix  = opts.YPadding;
boxW_pix  = opts.BoxWidth;
textX_pix = boxW_pix + xPad_pix;

% get parent figure
hParentFig = ancestor(hParentAx,'figure');

% save old units
oldParentAxUnits  = hParentAx.Units;
oldParentFigUnits = hParentFig.Units;

% for parent axes and figure set units and get normalized position
hParentAx.Units  = 'normalized';
hParentFig.Units = 'pixels';

%% Legend scaling: width and height

% create new axis under the created panel
axLeg = axes('Parent',hParentFig,'Position',[0.5 0.5 0.2 0.2],'YDir','reverse',...
    'Visible','on','Box','on','XLim',[0 1],'YLim',[0 1],...
    'XTick',[],'YTick',[]);

[xRange_pix, yRange_pix, rowH_pix, boxH_pix, ySpacer_pix, MaxTextW_pix_all, MaxTextH_pix] = findLegendSizePix(axLeg,entries,opts);
fixedWH_pix = [xRange_pix, yRange_pix];

% -------- set legend width and height in pixels ---------------
% set width and height in pixels
axLeg.Units = 'pixels';
axLeg.Position(3:4) = fixedWH_pix;

%% ==================== scale size relative to figure =====================
if opts.RelScaling
    % use axLeg.Units = pixels here to prevent rescaling of legend upon
    % rescaling of figure
    
    % get position of legend in pixels
    pos = getPos(axLeg,'pixels');

    % calculate relative width of legend
    relW  = pos(3)/getPos(hParentAx,'pixels',3);
    if relW > opts.MaxRelWidth
        % update width of figure - initial step
        hParentFig.Position(3) = hParentFig.Position(3)*relW/opts.MaxRelWidth;
        drawnow;
        
        % update width of figure - incremental changes
        relW  = pos(3)/getPos(hParentAx,'pixels',3);
        cW = 1; % counter
        while abs(relW/opts.MaxRelWidth-1) >= 0.01 && cW <= 100
            if cW < 3
                hParentFig.Position(3) = hParentFig.Position(3)*relW/opts.MaxRelWidth;
            elseif relW > opts.MaxRelWidth
                hParentFig.Position(3) = hParentFig.Position(3)+1;
            else
                hParentFig.Position(3) = hParentFig.Position(3)-1;
            end
            drawnow limitrate;
            relW  = pos(3)/getPos(hParentAx,'pixels',3);
            cW = cW + 1;
        end
    end

    % calculate relative height of le
    relH = pos(4)/getPos(hParentAx,'pixels',4);
    if relH > opts.MaxRelHeight
        % update height of figure - initial step
        hParentFig.Position(4) = hParentFig.Position(4)*relH/opts.MaxRelHeight;
        drawnow;

        % update height of figure - incremental changes
        relH  = pos(4)/getPos(hParentAx,'pixels',4);
        cH = 1; % counter
        while (abs(relH/opts.MaxRelHeight-1) >= 0.01 || relH > 1) && cH <= 100
            if cH < 3
                hParentFig.Position(4) = hParentFig.Position(4)*relH/opts.MaxRelHeight;
            elseif relH > opts.MaxRelHeight
                hParentFig.Position(4) = hParentFig.Position(4)+1;
            else
                hParentFig.Position(4) = hParentFig.Position(4)-1;
            end
            drawnow limitrate;
            relH  = pos(4)/getPos(hParentAx,'pixels',4);
        end
    end

end

%% ======================== Legend positioning ============================
% get position of legend in normalized coordinates

if isempty(pos0)
    % If no explicit Position, fall back to Location
    pos = CompassPositionHandler(axLeg,hParentAx,opts.Location);
else
    pos = axPosRel2figPosRel(hParentAx, pos0, 1:4);
end

% prevent making figure too small again
pos_init = getPos(axLeg,'normalized');
pos(3) = max(pos(3),pos_init(3));
pos(4) = max(pos(4),pos_init(4));

% set axLeg position
axLeg.Units = 'normalized';
axLeg.Position = pos;

%% Drawing entries

% Determine corresponding 'data' dimensions (i.e. in legend coordinates)
[MaxTextW,MaxTextH] = findMaxLegendTextWH(axLeg,entries,fontSize,'data');
MaxTextW_pix = max(MaxTextW_pix_all);
xPix2Data = MaxTextW / MaxTextW_pix;
yPix2Data = MaxTextH / MaxTextH_pix;

MaxTextW_all = MaxTextW_pix_all * xPix2Data;
rowH    = rowH_pix    * yPix2Data;
boxH    = boxH_pix    * yPix2Data;
ySpacer = ySpacer_pix * yPix2Data;
yPad    = yPad_pix    * yPix2Data;
xPad    = xPad_pix    * xPix2Data;
boxW    = boxW_pix    * xPix2Data;
textX   = textX_pix   * xPix2Data;

nEntries = numel(entries);
[i_all,j_all] = ind2sub([opts.rows,opts.cols],1:nEntries);
for k = 1:nEntries
    i = i_all(k);
    j = j_all(k);

    yCenter = (i-0.5)*rowH + (i-1)*ySpacer + yPad;
    y0 = yCenter - boxH/2;
    y1 = yCenter + boxH/2;

    xBase = sum(MaxTextW_all(1:j-1)) + (textX+2*xPad)*(j-1);

    % Patch
    patch(axLeg, [0 boxW boxW 0]+xPad + xBase, [y0 y0 y1 y1], entries(k).Color, ...
        'FaceAlpha', entries(k).Alpha, 'EdgeColor','none');
    
    % Separate the line style from marker
    [lineStyle, marker] = separateLineMarker(entries(k).LineStyle);

    % Line
    line(axLeg, [0 boxW]+xPad + xBase, [yCenter yCenter], ...
        'Color', entries(k).Color, ...
        'LineStyle', lineStyle, ...
        'LineWidth', entries(k).LineWidth);
    if ~isempty(marker)
        line(axLeg,boxW/2+xPad + xBase, yCenter, ...
        'Color', entries(k).Color, ...
        'Marker',marker,...
        'LineWidth', entries(k).LineWidth);
    end

    % Label
    text(axLeg, xPad+textX + xBase, yCenter, entries(k).Text, ...
        'Interpreter','latex', ...
        'VerticalAlignment','middle', ...
        'FontSize',fontSize);
end

%% --- Move legend to stay inside parent axes plotting box ---

% get the max/min relative position of the legend
set(hParentAx,'Units','normalized');
innerBox = get(hParentAx,'Position');

% legend position
pos = get(axLeg,'Position');

% Clamp legend inside parent axes/figure box
pos(1) = max(innerBox(1), min(pos(1), innerBox(1)+innerBox(3)-pos(3)));
pos(2) = max(innerBox(2), min(pos(2), innerBox(2)+innerBox(4)-pos(4)));
axLeg.Position = pos;

%% handling resizing
% set units to pixel (fixes width & height)
set(axLeg,'Units', 'pixels');

%% Reset units
hParentAx.Units  = oldParentAxUnits;
hParentFig.Units = oldParentFigUnits;

end

%% Helper functions
function posFig = CompassPositionHandler(hLeg,hAx,Location)

% get positions in normalized coordinates
posLeg = getPos(hLeg,'normalized');
posAx  = getPos(hAx, 'normalized');

% get width & height of legend w.r.t. axis
w = posLeg(3)/posAx(3);
h = posLeg(4)/posAx(4);

% specify position of legend w.r.t. axis
switch Location
    case "northwest"
        posLegInAx_xy = [0.05 0.95-h];
    case "northeast"
        posLegInAx_xy = [0.95-w 0.95-h];
    case "southwest"
        posLegInAx_xy = [0.05 0.05];
    case "southeast"
        posLegInAx_xy = [0.95-w 0.05];
    case "north"
        posLegInAx_xy = [0.5-w/2 0.95-h];
    case "south"
        posLegInAx_xy = [0.5-w/2 0.05];
    case "east"
        posLegInAx_xy = [0.95-w 0.5-h/2];
    case "west"
        posLegInAx_xy = [0.05 0.5-h/2];
    case "center"
        posLegInAx_xy = [0.5-w/2 0.5-h/2];
end
posLegInAx = [posLegInAx_xy w h];

% convert to position of legend w.r.t. figure
posFig = axPosRel2figPosRel(hAx, posLegInAx, 1:2);
posFig = [posFig posLeg(3) posLeg(4)];

end

function [MaxTextW,MaxTextH] = findMaxLegendTextWH(axLeg,entries,fontSize,Units)
htmp = text(axLeg, 0, 0, '', 'Interpreter','latex', ...
    'FontSize',fontSize,'Visible','off','Units',Units);
MaxTextW = 0;
MaxTextH = 0;
for k = 1:numel(entries)
    htmp.String = entries(k).Text;
    drawnow limitrate % faster than plain drawnow
    ext = htmp.Extent;
    MaxTextW = max(MaxTextW, ext(3));
    MaxTextH = max(MaxTextH, ext(4));
end
delete(htmp);
end

function posFig = axPosRel2figPosRel(hAx, posLeg, idx)
%AXPOS2FIGPOS Convert normalized position in axes to normalized position in figure
%
% posFig = axPos2figPos(axHandle, posAx)
%
% Inputs:
%   hAx      : handle to the reference axes
%   posLeg   : [x y w h] in normalized units w.r.t. the axes
%   idx      : specifies which of x, y, w, h to provide as output
%
% Output:
%   posFig   : [x y w h] in normalized units w.r.t. the parent axes

arguments
    hAx (1,1) matlab.graphics.axis.Axes
    posLeg (1,4) double
    idx double = [];
end
if ~isempty(idx)
    validateattributes(idx,{'numeric'},{'vector','real','finite','>=',1,'<=',4})
    if numel(idx) > 4
        error('idx may have only 4 indices to specify either, x, y, w or h');
    end
    idx = reshape(idx,1,numel(idx));
end

% Work in normalized units
posAx = getPos(hAx,'normalized'); % [x y w h] of axes in *figure normalized* coords

% Convert [x,y] from axes-normalized ? figure-normalized
xFig = posAx(1) + posLeg(1) * posAx(3);
yFig = posAx(2) + posLeg(2) * posAx(4);

% Convert [w,h] from axes-normalized ? figure-normalized
wFig = posLeg(3) * posAx(3);
hFig = posLeg(4) * posAx(4);

% Keep width/height unchanged
posFig = [xFig, yFig, wFig, hFig];

% output only a selection
posFig = posFig(idx);

end

function [xRange_pix, yRange_pix,rowH_pix,boxH_pix,ySpacer_pix,MaxTextW_pix_all,MaxTextH_pix] = findLegendSizePix(axLeg,entries,opts)
rowHfac   = opts.RowHeightFactor;
boxHfac   = opts.BoxHeightFactor;
ySpaceFac = opts.YSpacingFactor;
fontSize  = opts.FontSize;
xPad_pix  = opts.XPadding;
yPad_pix  = opts.YPadding;
boxW_pix  = opts.BoxWidth;
textX_pix = boxW_pix + xPad_pix;


% measure text height and width (in pixels) for each column in the legend
nEntries = numel(entries);
[~,j] = ind2sub([opts.rows,opts.cols],1:nEntries);
MaxTextW_pix_all = nan(1,opts.cols);
MaxTextH_pix_all = nan(1,opts.cols);
for kc = 1:opts.cols
    [MaxTextW_pix_all(kc),MaxTextH_pix_all(kc)] = findMaxLegendTextWH(axLeg,entries(j==kc),fontSize,'pixels');
end
MaxTextH_pix = max(MaxTextH_pix_all);

% way to calculate the width
xRangeFun = @(textX,MaxTextWidth_sum,xPad) textX*opts.cols + MaxTextWidth_sum + (opts.cols-1)*2*xPad + 2*xPad;

% way to calculate the height
yRangeFun = @(rowH,ySpacer,yPad) opts.rows*rowH + (opts.rows-1)*ySpacer +2*yPad;

% -> calculate dependent constants
rowH_pix    = rowHfac   * MaxTextH_pix;
boxH_pix    = boxHfac   * MaxTextH_pix;
ySpacer_pix = ySpaceFac * MaxTextH_pix;

% determine desired legend width and height in pixels
xRange_pix = xRangeFun(textX_pix,sum(MaxTextW_pix_all),xPad_pix);
yRange_pix = yRangeFun(rowH_pix,ySpacer_pix,yPad_pix);
end

% get position in described units but don't change units thereafter
function Pos = getPos(hobj,Units,idx)
if nargin < 3
    idx = 1:4;
end
validateattributes(idx, {'numeric'},{'vector','real','finite','>=',1,'<=',4});

oldUnits = hobj.Units;
hobj.Units = Units;

Pos = hobj.Position;
% if isa(hobj,"matlab.graphics.axis.Axes")
%     if isa(get(hobj, 'parent'),"matlab.graphics.layout.TiledChartLayout")
%         % % disable warning for axes on TiledChartLayout
%         % warning_id = 'MATLAB:handle_graphics:Layout:NoPositionSetInTiledChartLayout';
%         % warning('off',warning_id);
%         Pos = hobj.Position;
%     else
%         Pos = plotboxpos(hobj); % would require this function from the File Exchange
%     end    
% else
%     Pos = hobj.Position;
% end

Pos = Pos(idx);
hobj.Units = oldUnits;
end

function opts = validateOptions(opts,entries)
% calculate number of rows
opts.rows = ceil(numel(entries)/opts.cols);

% Validate Position if the user provided one (empty means "not provided")
if ~isempty(opts.Position)
    validateattributes(opts.Position, {'numeric'}, ...
        {'vector','numel',4,'real','finite','>=',0,'<=',1}, ...
        mfilename, 'Position');
    opts.Position = reshape(opts.Position, 1, 4);
    if ~isempty(opts.MaxRelWidth) && opts.Position(3) > opts.MaxRelWidth
        error('Position(3) = %f, but must be less than or equal to the MaxRelWidth parameter (=%f)', ...
            opts.Position(3), opts.MaxRelWidth)
    elseif isempty(opts.MaxRelWidth)
        opts.MaxRelWidth = opts.Position(3);
    end
    if ~isempty(opts.MaxRelHeight) && opts.Position(4) > opts.MaxRelHeight
        error('Position(4) = %f, but must be less than or equal to the MaxRelHeight parameter (=%f)', ...
            opts.Position(4), opts.MaxRelHeight)
    elseif isempty(opts.MaxRelHeight)
        opts.MaxRelHeight = opts.Position(4);
    end
end

% set defaults for MaxRelWidth & MaxRelHeight
if isempty(opts.MaxRelWidth)
    opts.MaxRelWidth = 0.2 * opts.cols;
end
if isempty(opts.MaxRelHeight)
    opts.MaxRelHeight = 0.9;
end

% Validate relative legend dimensons
if ~isempty(opts.MaxRelWidth)
    validateattributes(opts.MaxRelWidth, {'numeric'}, ...
        {'scalar','real','finite','>=',0,'<=',0.9},mfilename,'MaxRelWidth');
end
if ~isempty(opts.MaxRelWidth)
    validateattributes(opts.MaxRelHeight, {'numeric'}, ...
        {'scalar','real','finite','>=',0,'<=',0.9},mfilename,'MaxRelHeight');
end

end

function [lineStyle, marker] = separateLineMarker(styleStr)
    % separateLineMarker - Separates line style from marker type
    %
    % Syntax: [lineStyle, marker] = separateLineMarker(styleStr)
    %
    % Inputs:
    %   styleStr - String containing line style and/or marker (e.g., '-o', '--x', ':')
    %
    % Outputs:
    %   lineStyle - Line style ('-', '--', ':', '-.', or 'none')
    %   marker - Marker type ('o', 'x', '+', etc., or 'none')
    
    % Initialize outputs
    lineStyle = 'none';
    marker = 'none';
    
    if isempty(styleStr)
        return;
    end
    
    % Define valid line styles (order matters - check longer ones first)
    lineStyles = {'--', '-.', '-', ':'};
    
    % Define valid markers
    markers = {'+', 'o', '*', '.', 'x', 's', 'd', '^', 'v', '>', '<', ...
               'p', 'h', 'square', 'diamond', 'pentagram', 'hexagram'};
    
    % Check for line style
    for i = 1:length(lineStyles)
        if contains(styleStr, lineStyles{i})
            lineStyle = lineStyles{i};
            % Remove line style from string
            styleStr = strrep(styleStr, lineStyles{i}, '');
            break;
        end
    end
    
    % Check for marker in remaining string
    for i = 1:length(markers)
        if contains(styleStr, markers{i})
            marker = markers{i};
            break;
        end
    end
    
    % If nothing was found, default to solid line
    if strcmp(lineStyle, 'none') && strcmp(marker, 'none')
        lineStyle = '-';
    end
end