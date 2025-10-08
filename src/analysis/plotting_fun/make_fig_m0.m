function [fig1,ax1,axLeg1] = make_fig_m0(m0, pX, useFillCases, noPlotCases,opts)
arguments
    m0 struct
    pX (1,1) double {mustBeInteger,mustBePositive}
    useFillCases cell
    noPlotCases cell
    opts.CrameriColors (1,:) char = 'roma';  % possible colors: see the 'crameri' command
    opts.FigPos (1,4) double {mustBeReal,mustBeFinite} = [50 50 1000 600]; % figure position
    opts.fillAlpha (1,1) double {mustBePositive,mustBeReal,mustBeFinite,...
                                 mustBeLessThanOrEqual(opts.fillAlpha,1)} = 0.25
    opts.LineWidth (1,1) double {mustBeReal,mustBeFinite,mustBePositive} = 2
    opts.LegLocations (1,:) = "northeast"; % parsing done by mustBeValidLegLocation
    opts.LegCols (1,2) double {mustBeFinite,mustBeReal,mustBeGreaterThan(opts.LegCols,0),mustBeInteger} = [1, 1];
    opts.FS_Label (1,1) double {mustBeReal,mustBeFinite,mustBePositive} = 15;
    opts.FS_Legend (1,1) double {mustBeReal,mustBeFinite,mustBePositive} = 12;
end
%MAKE_FIG_M0  Prepare and plot IV trajectories with bounds for a given index.
%
%   make_fig_m0(m0, pX, opts)
%
%   Inputs:
%     m0   Struct with fields Uf and Yf, each containing IV data.
%     pX   Index for the iX string (integer).
%     useFillCases Cell array of case names for which to fill percentile bounds.
%     noPlotCases  Cell array of case names to exclude from plotting.
%     opts Options struct passed to plot_IV_trajectories.

cColors   = opts.CrameriColors;
Position  = opts.FigPos;
fillAlpha = opts.fillAlpha;
LineWidth = opts.LineWidth;
LegsLoc   = mustBeValidLegLocation(opts.LegLocations);
LegsCols  = opts.LegCols;
fontSize_xyLabel = opts.FS_Label;
fontSize_legend = opts.FS_Legend;

%% Preliminaries - iv names & colours
% get number of field names in Uf and Yf -> # of colors needed
Uf_ivs = sort(fieldnames(m0.Uf));
Yf_ivs = sort(fieldnames(m0.Yf));

% only plot cases not in noPlotCases
Uf_ivs = setdiff(Uf_ivs,noPlotCases);
Yf_ivs = setdiff(Yf_ivs,noPlotCases);

% moving cases around for better visibility
msk_uiv6 = contains(Uf_ivs,'iv6'); Uf_ivs = [Uf_ivs(msk_uiv6); Uf_ivs(~msk_uiv6)]; % move iv6 to beginning
msk_yiv6 = contains(Yf_ivs,'iv6'); Yf_ivs = [Yf_ivs(msk_yiv6); Yf_ivs(~msk_yiv6)];

ivName2ColorIdx = @(s) str2double(regexp(s,'\d+','match','once')); % extract number from iv name
Cidxs = unique(cellfun(@(s) ivName2ColorIdx(s), [Uf_ivs;Yf_ivs]));
Cidxs = setdiff(Cidxs,6); % remove case 6 since these are plotted in black
msk_1 = Cidxs == 1;
msk_2 = Cidxs == 2;
Cidxs = [Cidxs(msk_2); Cidxs(~(msk_1 | msk_2)); Cidxs(msk_1);];
% if Cidxs(1) == 1
%     Cidxs = circshift(Cidxs,-1); % clearly separate cases 1 & 2
% end
% Cidxs = flipud(Cidxs);
numColors = numel(Cidxs);

load("pdir.mat",'pdir');
addpath(fullfile(pdir,'bin','external','crameri_colours'))
cCram = crameri(cColors, numColors); % colors

% Define line styles to cycle through
lineStyles = {'-','--',':','-.'};

iX_str = sprintf('iX%d',pX);
% select the 25th and 75th percentile -> 6 & 16th of 21
% see size(m0.Uf.iv1.iX1.pctiles,3);

fig1 = figure("Position",Position);
tl1 = tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

%% ============================= Uf IVs =================================
ax1(1) = nexttile(tl1,1);
Uf_lb_fs = replace(Uf_ivs,'iv','l');
Uf_ub_fs = replace(Uf_ivs,'iv','u');
Uf_av_fs = replace(Uf_ivs,'iv','m');
num_Uf_ivs = numel(Uf_ivs);

u_iv = struct;
for kIV = 1:num_Uf_ivs
    lb_name = Uf_lb_fs{kIV};
    ub_name = Uf_ub_fs{kIV};
    av_name = Uf_av_fs{kIV};
    iv_name = Uf_ivs{kIV};
    
    mean_vals = m0.Uf.(iv_name).(iX_str).median;
    lb_vals = squeeze(m0.Uf.(iv_name).(iX_str).pctiles(:,:,6));
    ub_vals = squeeze(m0.Uf.(iv_name).(iX_str).pctiles(:,:,16));

    cIdx = ivName2ColorIdx(iv_name);
    cIdx2 = find(Cidxs == cIdx,1,'first');
    
    switch iv_name
        case 'iv1',  label = '1) open-loop IV';          col = cCram(cIdx2,:);
        case 'iv2a', label = '2a) opt. IV';              col = cCram(cIdx2,:);
        case 'iv2c', label = '2c) opt. IV + 2SLS';       col = cCram(cIdx2,:);
        case 'iv3c', label = '3c) LCF + 2SLS';           col = cCram(cIdx2,:);
        case 'iv4a', label = '4a) w/o Cz info';          col = cCram(cIdx2,:);
        case 'iv4c', label = '4c) w/o Cz info + 2SLS';   col = cCram(cIdx2,:);
        case 'iv5a', label = '5a) w/ Cz info';           col = cCram(cIdx2,:);
        case 'iv5c', label = '5c) w/ Cz info + 2SLS';    col = cCram(cIdx2,:);
        case 'iv6c', label = '6c) ref + 2SLS';           col = [0 0 0]; % black
        otherwise,  label = iv_name;                     col = [0 0 0];
    end
    switch iv_name
        case 'iv2a', thisStyle = '-';
        otherwise, thisStyle = lineStyles{mod(kIV-1,numel(lineStyles))+1}; % cycle
    end
    if any(ismember(iv_name,useFillCases))
        useFill = true;
    else
        useFill = false;
    end

    plotLineWithFill(mean_vals, lb_vals, ub_vals, label, Color=col, ...
        FaceAlpha = fillAlpha, LineStyle=thisStyle, LineWidth=LineWidth,...
        XScale='linear',YScale='linear',useFill = useFill);
    hold on;

    u_entries(kIV) = struct('Color',col, 'Alpha',fillAlpha*useFill, ...
                    'LineStyle',thisStyle, 'LineWidth',LineWidth, ...
                    'Text',label,'FontSize',fontSize_legend);
end
grid on;
ylabel('$u_k$','Interpreter','latex','FontSize',fontSize_xyLabel);
axLeg1(1) = customLegend(u_entries,ax1(1),cols=LegsCols(1),Location=LegsLoc(1),RelScaling=true);

%% ============================= Yf IVs ================================
ax1(2) = nexttile(tl1,2);

Yf_lb_fs = replace(Yf_ivs,'iv','l');
Yf_ub_fs = replace(Yf_ivs,'iv','u');
Yf_av_fs = replace(Yf_ivs,'iv','m');
num_Yf_ivs = numel(Yf_ivs);

y_iv = struct;
for kIV = 1:num_Yf_ivs
    lb_name = Yf_lb_fs{kIV};
    ub_name = Yf_ub_fs{kIV};
    av_name = Yf_av_fs{kIV};
    iv_name = Yf_ivs{kIV};
    
    mean_vals = m0.Yf.(iv_name).(iX_str).median;
    lb_vals = squeeze(m0.Yf.(iv_name).(iX_str).pctiles(:,:,6));
    ub_vals = squeeze(m0.Yf.(iv_name).(iX_str).pctiles(:,:,16));

    cIdx = ivName2ColorIdx(iv_name);
    cIdx2 = find(Cidxs == cIdx,1,'first');

    switch iv_name
        case 'iv2b', label = '2b) opt. IV';             col = cCram(cIdx2,:);
        case 'iv4b', label = '4b) w/o Cz info';         col = cCram(cIdx2,:);
        case 'iv5b', label = '5b) w/ Cz info';          col = cCram(cIdx2,:);
        case 'iv3a', label = '3a) LCF-IV ($\Theta$)';   col = cCram(cIdx2,:);
        case 'iv6a', label = '6a) ref';                 col = [0 0 0]; % black
        otherwise, label = iv_name;                     col = [0 0 0];
    end
    thisStyle = lineStyles{mod(kIV-1,numel(lineStyles))+1}; % cycle
    if any(ismember(iv_name,useFillCases))
        useFill = true;
    else
        useFill = false;
    end

    plotLineWithFill(mean_vals, lb_vals, ub_vals, label, Color=col, ...
        FaceAlpha = fillAlpha, LineStyle=thisStyle, LineWidth=LineWidth,...
        XScale='linear',YScale='linear',useFill = useFill);
    hold on;

    y_entries(kIV) = struct('Color',col, 'Alpha',fillAlpha*useFill, ...
                    'LineStyle',thisStyle, 'LineWidth',LineWidth, ...
                    'Text',label,'FontSize',fontSize_legend);
end
grid on;
ylabel('$y_k$','Interpreter','latex','FontSize',fontSize_xyLabel);
xlabel('Time step', 'Interpreter','latex','FontSize',fontSize_xyLabel);
axLeg1(2) = customLegend(y_entries,ax1(2),cols=LegsCols(2),Location=LegsLoc(2),RelScaling=true);

linkaxes(ax1,'x');

end

%% Helper functions
function val = mustBeValidLegLocation(val)
% mustBeValidLegLocation  Validate legend location(s) with replication rule
%
%   Rules:
%   - Accepts string arrays, char, or cell arrays containing only char/string.
%   - Must have 1 or 2 elements.
%   - If only 1 element is provided, it is replicated to 2 elements.
%   - If input is a cell array, it is converted to a string array.
%   - Each element must be a valid legend location.

    % Allowed values
    validLocations = ["northeast","northwest","southeast","southwest", ...
                      "east","west","north","south"];

    % Check number of elements
    n = numel(val);
    if n < 1 || n > 2
        error("Leg.Location must have 1 or 2 elements, but has %d.", n);
    end

    % Normalize input type
    if ischar(val)
        val = string(val);

    elseif iscell(val)
        % Check that every element is char or string
        if ~all(cellfun(@(x) ischar(x) || isstring(x), val))
            error("Leg.Location cell array must only contain strings or character vectors.");
        end

        % Convert all chars to strings
        val = cellfun(@string, val);

        % Convert cell -> string array
        val = string(val);

    elseif ~isstring(val)
        error("Leg.Location must be a string array, char, or cell array of char/string.");
    end

    % Validate each element
    for i = 1:n
        if ~any(validLocations == val(i))
            error("Invalid legend location: '%s'. Allowed values are: %s", ...
                val(i), strjoin(validLocations,", "));
        end
    end

    % Replicate if single element
    if n == 1
        val = [val, val];
    end
end