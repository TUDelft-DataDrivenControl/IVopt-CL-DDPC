function axLeg = customLegend(entries, parentFig, opts)
%CUSTOMLEGEND Create a custom legend with patch+line+LaTeX text
% The legend auto-resizes/repositions itself to fit all text inside.
%
%   axLeg = customLegend(parentFig, entries)
%   axLeg = customLegend(parentFig, entries, Name,Value,...)
%
% Inputs:
%   parentFig : figure or axes handle
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

arguments
    entries (1,:) struct
    parentFig {mustBeA(parentFig,["matlab.ui.Figure","matlab.graphics.axis.Axes"])} = gca
    opts.Position = []; % explicit [x y w h], empty = auto
    opts.XPadding (1,1) double = 5  % [pixels]
    opts.YPadding (1,1) double = 5  % [pixels]
    opts.BoxWidth (1,1) double = 30 % [pixels]
    opts.RowHeightFactor (1,1) double = 1
    opts.BoxHeightFactor (1,1) double = 1
    opts.YSpacingFactor (1,1) double = 0.2
    opts.FontSize (1,1) double = 12
    opts.Location (1,1) string {mustBeMember(opts.Location, ...
        ["northwest","northeast","southwest","southeast","north","south","east","west","center"])} = "northeast"
end
nEntries = numel(entries); % number of legend entries

yRangeFun = @(rowH,ySpacer,yPad) nEntries*rowH + (nEntries-1)*ySpacer +2*yPad;
xRangeFun = @(textX,MaxTextWidth,xPad) textX + MaxTextWidth + 2*xPad;
% Validate Position if the user provided one (empty means "not provided")
if ~isempty(opts.Position)
    validateattributes(opts.Position, {'numeric'}, ...
        {'vector','numel',4,'real','finite','>=',0,'<=',1}, ...
        mfilename, 'Position');
    opts.Position = reshape(opts.Position, 1, 4);
end
pos0      = opts.Position;
rowHfac   = opts.RowHeightFactor;
boxHfac   = opts.BoxHeightFactor;
ySpaceFac = opts.YSpacingFactor;
fontSize  = opts.FontSize;
xPad_pix  = opts.XPadding;
yPad_pix  = opts.YPadding;
boxW_pix  = opts.BoxWidth;
textX_pix = boxW_pix + xPad_pix;

% --- Resolve legend position ---
if isempty(pos0)
    % If no explicit Position, fall back to Location
    pos = CompassPositionHandler(opts.Location);    
else
    pos = pos0;
end

% If parent is an axes, get its figure
if isa(parentFig,'matlab.graphics.axis.Axes')
    parentAx = parentFig;
    parentIsAx = true;
    parentFig = ancestor(parentFig,'figure');
else
    parentIsAx = false;
end

%% Legend scaling: width and height
axLeg = axes('Parent',parentFig,'Position',pos,'YDir','reverse',...
    'Visible','on','Box','on','XLim',[0 1],'YLim',[0 1],...
    'XTick',[],'YTick',[]);

% measure text height and width (in pixels)
[MaxTextW_pix,MaxTextH_pix] = findMaxTextWH(axLeg,entries,fontSize,'pixels');

% -> calculate dependent constants
rowH_pix    = rowHfac   * MaxTextH_pix;
boxH_pix    = boxHfac   * MaxTextH_pix;
ySpacer_pix = ySpaceFac * MaxTextH_pix;

% determine desired legend width and height in pixels
xRange_pix = xRangeFun(textX_pix,MaxTextW_pix,xPad_pix);
yRange_pix = yRangeFun(rowH_pix,ySpacer_pix,yPad_pix);

% set legend width and height in pixels
axLeg.Units = 'pixels';
axLeg.Position(3:4) = [xRange_pix, yRange_pix];
axLeg.Units = 'normalized';
pos = axLeg.Position;
if isempty(pos0)
    % If no explicit Position, fall back to Location
    pos = CompassPositionHandler(opts.Location,pos);
end

% set position relative to Axis or Figure
if parentIsAx
    pos = axPos2figPos(parentAx, pos);
end
axLeg.Position = pos;

% Determine corresponding 'data' dimensions (i.e. in legend coordinates)
[MaxTextW,MaxTextH] = findMaxTextWH(axLeg,entries,fontSize,'data');
xPix2Data = MaxTextW / MaxTextW_pix;
yPix2Data = MaxTextH / MaxTextH_pix;

rowH    = rowH_pix    * yPix2Data;
boxH    = boxH_pix    * yPix2Data;
ySpacer = ySpacer_pix * yPix2Data;
yPad    = yPad_pix    * yPix2Data;
xPad    = xPad_pix    * xPix2Data;
boxW    = boxW_pix    * xPix2Data;
textX   = textX_pix   * xPix2Data;

%% Drawing entries
for k = 1:numel(entries)
    yCenter = (k-0.5)*rowH + (k-1)*ySpacer + yPad;
    y0 = yCenter - boxH/2;
    y1 = yCenter + boxH/2;

    % Patch
    patch(axLeg, [0 boxW boxW 0]+xPad, [y0 y0 y1 y1], entries(k).Color, ...
        'FaceAlpha', entries(k).Alpha, 'EdgeColor','none');

    % Line
    line(axLeg, [0 boxW]+xPad, [yCenter yCenter], ...
        'Color', entries(k).Color, ...
        'LineStyle', entries(k).LineStyle, ...
        'LineWidth', entries(k).LineWidth);

    % Label
    text(axLeg, xPad+textX, yCenter, entries(k).Text, ...
        'Interpreter','latex', ...
        'VerticalAlignment','middle', ...
        'FontSize',12);
end

%% --- Move legend to stay inside parent axes plotting box ---

% get the max/min relative position of the legend
if parentIsAx
    pAxUnits = parentAx.Units;
    set(parentAx,'Units','normalized');
    innerBox = get(parentAx,'Position');
else
    innerBox = [0 0 1 1];
end

% legend position
pos = get(axLeg,'Position');

% Clamp legend inside parent axes/figure box
pos(1) = max(innerBox(1), min(pos(1), innerBox(1)+innerBox(3)-pos(3)));
pos(2) = max(innerBox(2), min(pos(2), innerBox(2)+innerBox(4)-pos(4)));
axLeg.Position = pos;

% reset units
if parentIsAx
    set(parentAx,'Units',pAxUnits);
end

end

%% Helper functions
function pos = CompassPositionHandler(Location,pos)

if nargin < 2
    % base width/height (could be tuned)
    w = 0.25;
    h = 0.2;
else
    w = pos(3);
    h = pos(4);
end

switch Location
    case "northwest"
        pos = [0.05 0.95-h w h];
    case "northeast"
        pos = [0.95-w 0.95-h w h];
    case "southwest"
        pos = [0.05 0.05 w h];
    case "southeast"
        pos = [0.95-w 0.05 w h];
    case "north"
        pos = [0.5-w/2 0.95-h w h];
    case "south"
        pos = [0.5-w/2 0.05 w h];
    case "east"
        pos = [0.95-w 0.5-h/2 w h];
    case "west"
        pos = [0.05 0.5-h/2 w h];
    case "center"
        pos = [0.5-w/2 0.5-h/2 w h];
end
end

function [MaxTextW,MaxTextH] = findMaxTextWH(axLeg,entries,fontSize,Units)
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

function posFig = axPos2figPos(axHandle, posAx)
%AXPOS2FIGPOS Convert normalized position in axes to normalized position in figure
%
% posFig = axPos2figPos(axHandle, posAx)
%
% Inputs:
%   axHandle : handle to the reference axes
%   posAx    : [x y w h] in normalized units w.r.t. the axes
%
% Output:
%   posFig   : [x y w h] in normalized units w.r.t. the parent figure
%
% Note: width/height (w,h) are passed through unchanged.

    arguments
        axHandle (1,1) matlab.graphics.axis.Axes
        posAx (1,4) double {mustBeGreaterThanOrEqual(posAx,0),mustBeLessThanOrEqual(posAx,1)}
    end

    % Store old units
    oldAxUnits = axHandle.Units;
    oldFigUnits = axHandle.Parent.Units;

    % Work in normalized units
    axHandle.Units = 'normalized';
    axPosFig = axHandle.Position; % [x y w h] of axes in *figure normalized* coords

    % Convert [x,y] from axes-normalized → figure-normalized
    xFig = axPosFig(1) + posAx(1) * axPosFig(3);
    yFig = axPosFig(2) + posAx(2) * axPosFig(4);

    % Keep width/height unchanged
    posFig = [xFig, yFig, posAx(3), posAx(4)];

    % Restore units
    axHandle.Units = oldAxUnits;
    axHandle.Parent.Units = oldFigUnits;
end