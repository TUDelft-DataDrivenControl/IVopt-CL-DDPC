function h = plotLineWithFill(y, yLower, yUpper, description, varargin)
% plotLineWithFill plots a line with shaded bounds and a single legend
% entry that shows both the line and the patch together.
%
% Usage:
%   plotLineWithFill(y, yLower, yUpper, description)
%   plotLineWithFill(y, yLower, yUpper, description, 'x', xVals)
%   plotLineWithFill(y, yLower, yUpper, description, 'LineStyle','--','FaceAlpha',0.5)
%   plotLineWithFill(y, yLower, yUpper, description, 'XScale','log')
%
% Inputs:
%   y          - main line values
%   yLower     - lower bound (same length as y)
%   yUpper     - upper bound (same length as y)
%   description- string used for legend
%
% Optional name-value pairs:
%   'x'        - x values (default = 1:length(y))
%   'LineStyle' - line style for main line (default = '-')
%   'FaceAlpha' - transparency of fill (default = 0.25)
%   'Color'     - RGB triplet or color char (default = next color in ColorOrder)
%   'LineWidth' - line width (default = 1.5)
%   'XScale'    - 'linear' (default) or 'log'
%   'YScale'    - 'linear' (default) or 'log'
%   'useFill'   - true (default) or false
%
% Output:
%   h - structure with handles to line and fill

    % Ensure column-to-row formatting
    y      = y(:).';
    yLower = yLower(:).';
    yUpper = yUpper(:).';

    % Parse inputs
    p = inputParser;
    addParameter(p, 'LineStyle', '-', @ischar);
    addParameter(p, 'FaceAlpha', 0.25, @(a)isnumeric(a)&&isscalar(a));
    addParameter(p, 'Color', [], @(c)(ischar(c) || (isnumeric(c)&&numel(c)==3)));
    addParameter(p, 'LineWidth', 1.5, @(a)isnumeric(a)&&isscalar(a));
    addParameter(p, 'x', 1:length(y), @isnumeric);
    addParameter(p, 'XScale', 'linear', @(s)ischar(s) && any(strcmpi(s,{'linear','log'})));
    addParameter(p, 'YScale', 'linear', @(s)ischar(s) && any(strcmpi(s,{'linear','log'})));
    addParameter(p, 'useFill', true, @islogical);
    parse(p, varargin{:});

    % Extract results
    ls     = p.Results.LineStyle;
    fa     = p.Results.FaceAlpha;
    c      = p.Results.Color;
    lw     = p.Results.LineWidth;
    x      = p.Results.x(:).';
    xScale = lower(p.Results.XScale);
    yScale = lower(p.Results.YScale);
    useFill= p.Results.useFill;

    % Default x if empty
    if isempty(x)
        x = 1:length(y);
    end

    % Get axes & color order
    ax = gca;
    if isempty(c)
        co = get(ax,'ColorOrder');
        nLines = numel(findall(ax,'type','line'));
        c = co(mod(nLines, size(co,1))+1, :);
    end

    % --- plot patch and line ---
    if useFill
        h.fill = fill(ax,[x fliplr(x)], [yLower fliplr(yUpper)], c, ...
                      'FaceAlpha', fa, 'EdgeColor','none', ...
                      'HandleVisibility','off');%'DisplayName', description);
        hold on;
    end
    h.line = plot(ax, x, y, 'LineStyle', ls, 'Color', c, 'LineWidth', lw, ...
                  'HandleVisibility','off'); % hide line in legend

    % --- set axis scaling ---
    set(ax, 'XScale', xScale);
    set(ax, 'YScale', yScale);
end