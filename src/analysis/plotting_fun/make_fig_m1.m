function [fig2,ax2,axLeg2] =make_fig_m1(m1,data_type,X_all,opts)
arguments
    m1 struct
    data_type (1,:) char {mustBeMember(data_type,{'N','Re','p'})}
    X_all double {mustBeVector}
    opts.N (1,1) double {mustBeFinite,mustBeReal,mustBeInteger,mustBePositive};
    opts.fontSize (1,1) double = 15;
    opts.CrameriColors char = 'roma';  % possible colors: see the 'crameri' command
    opts.FigPos (1,4) double {mustBeReal,mustBeFinite} = [50 50 1000 600]; % figure position
    opts.Units = 'pixels'
    opts.fillAlpha (1,1) double {mustBeGreaterThanOrEqual(opts.fillAlpha,0),...
                                 mustBeLessThanOrEqual(opts.fillAlpha,1)} = 0.25
    opts.XScale char {mustBeMember(opts.XScale,["log","linear"])} = 'log';
    opts.YScale char {mustBeMember(opts.YScale,["log","linear"])} = 'linear';
    opts.LineWidth (1,1) double {mustBeFinite,mustBeReal,mustBePositive} = 2;
    opts.LegLocations (1,:) = "northeast"; % parsing done by mustBeValidLegLocation
    opts.LegCols (1,2) double {mustBeFinite,mustBeReal,mustBePositive,mustBeInteger} = [1, 1];
    opts.FS_Legend (1,1) double {mustBeReal,mustBeFinite,mustBePositive} = 12;
end

fontSize  = opts.fontSize;      % font size of x & y labels (scaled later for y due to orientation)
cColors   = opts.CrameriColors;
Position  = opts.FigPos;
fillAlpha = opts.fillAlpha;
XScale    = opts.XScale;
YScale    = opts.YScale;
LineWidth = opts.LineWidth;
LegsLoc   = mustBeValidLegLocation(opts.LegLocations);
LegsCols  = opts.LegCols;

% Determine x-axis
switch data_type
    case 'N'
        N_all = X_all;
        xlab = '$N$';
    case {'Re','p'}
        if ~isfield(opts,'N')
            error('Specify N used')
        end
        N_all = repmat(opts.N,1,numel(X_all));

        if strcmp(data_type,'Re')
            xlab = '$\Sigma_e$';
        else
            xlab = '$p$';
        end        
end

% --------------------- set colours and make figure -----------------------
try
    cCram = crameri(cColors, 9); % colors
catch ME
    msg = ME.message;
    msg = sprintf('Crameri color palette not recognized, resulting in:\n%s',msg);
    ME.message = msg;
   rethrow(ME)
end


set_xlabel_2   = @() xlabel(xlab,'interpreter','latex','FontSize',fontSize);
set_ylabel_2_1 = @() ylabel('$\frac{\|\Delta_j \tilde{U}_{\mathrm{f}}\|_\mathrm{F}}{\sqrt{N}}$',...
    'Interpreter','latex','Rotation',0,'FontSize',fontSize*1.25);
set_ylabel_2_2 = @() ylabel('$\frac{\|\Delta_j \tilde{Y}_{\mathrm{f}}\|_\mathrm{F}}{\sqrt{N}}$',...
    'Interpreter','latex','Rotation',0,'FontSize',fontSize*1.25);

fig2 = figure('Units',opts.Units,'Position',Position);
tl2 = tiledlayout(2,1,"TileSpacing",'tight','Padding','tight');
fig2.Units = 'pixels';
ax2(1) = nexttile();

% --------------------- plotting for Uf_ivs -------------------------------
% plotting Uf frobenius norm errors
% select the 25th and 75th percentile -> 6 & 16th of 21
% see size(m1.Uf.iv1.pctiles,2);

Uf_ivs     = fieldnames(m1.Uf);
Uf_ivs = setdiff(Uf_ivs,{'iv2a'}); % exclude case that is logically zero
Uf_lb_fs = replace(Uf_ivs,'iv','l');
Uf_ub_fs = replace(Uf_ivs,'iv','u');
Uf_av_fs = replace(Uf_ivs,'iv','m');
num_Uf_ivs = numel(Uf_ivs);

sqrt_Nall = diag(sqrt(N_all));
u_iv = struct;
for kIV = 1:num_Uf_ivs
    lb_name = Uf_lb_fs{kIV};
    ub_name = Uf_ub_fs{kIV};
    av_name = Uf_av_fs{kIV};
    iv_name = Uf_ivs{kIV};
    
    u_iv.(av_name) = sqrt_Nall\m1.Uf.(iv_name).median;
    u_iv.(lb_name) = sqrt_Nall\squeeze(m1.Uf.(iv_name).pctiles(:,6));
    u_iv.(ub_name) = sqrt_Nall\squeeze(m1.Uf.(iv_name).pctiles(:,16));
    
    % determine color
    switch iv_name
        case 'iv6c', label = '$\hat{\tilde{U}}_{\mathrm{f,6c}}$'; col = [0 0 0];       LineStyle = '-';
        case 'iv2c', label = '$\hat{\tilde{U}}_{\mathrm{f,2c}}$'; col = cCram(2,:);    LineStyle = '-.';
        case 'iv3c', label = '$\hat{\tilde{U}}_{\mathrm{f,3c}}$'; col = cCram(3,:);    LineStyle = '-';
        case 'iv4a', label = '$\hat{\tilde{U}}_{\mathrm{f,4a}}$'; col = cCram(4,:);    LineStyle = '-*';
        case 'iv4c', label = '$\hat{\tilde{U}}_{\mathrm{f,4c}}$'; col = cCram(5,:);    LineStyle = '-^';
        case 'iv5a', label = '$\hat{\tilde{U}}_{\mathrm{f,5a}}$'; col = cCram(6,:);    LineStyle = '-.*';
        case 'iv5c', label = '$\hat{\tilde{U}}_{\mathrm{f,5c}}$'; col = cCram(7,:);    LineStyle = '-.^';
        case 'iv1',  label = '$U_{\mathrm{f}}$';                  col = cCram(9,:);    LineStyle = '-';
        otherwise,   label = iv_name;                             col = [0 0 0];       LineStyle = '-';
    end

    plotLineWithFill(u_iv.(av_name),u_iv.(lb_name),u_iv.(ub_name), label, ...
        Color=col, FaceAlpha = fillAlpha, LineStyle=LineStyle, LineWidth=LineWidth,...
         x=X_all, XScale=XScale,YScale=YScale);
    hold on;

    u_entries(kIV) = struct('Color',col, 'Alpha',fillAlpha, ...
                    'LineStyle',LineStyle, 'LineWidth',LineWidth, ...
                    'Text',label);
end
set_ylabel_2_1();
grid on;
axLeg2(1) = customLegend(u_entries,ax2(1),Location=LegsLoc(1),cols=LegsCols(1),FontSize=opts.FS_Legend);

% --------------------- plotting for Yf_ivs -------------------------------
ax2(2) = nexttile(tl2,2);

Yf_ivs     = fieldnames(m1.Yf);
Yf_ivs = setdiff(Yf_ivs,{'iv2b'}); % exclude case that is logically zero
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
    
    y_iv.(av_name) = sqrt_Nall\m1.Yf.(iv_name).median;
    y_iv.(lb_name) = sqrt_Nall\squeeze(m1.Yf.(iv_name).pctiles(:,6));
    y_iv.(ub_name) = sqrt_Nall\squeeze(m1.Yf.(iv_name).pctiles(:,16));

    % determine color
    switch iv_name
        case 'iv4b', label = '$\hat{\tilde{Y}}_{\mathrm{f,4b}}$'; col = cCram(4,:);  LineStyle = '-';
        case 'iv5b', label = '$\hat{\tilde{Y}}_{\mathrm{f,5b}}$'; col = cCram(6,:);  LineStyle = '-.d';
        case 'iv3a', label = '$\Xi_f$' ;                          col = cCram(8,:);  LineStyle = '-x';
        case 'iv6a', label = '$W_f$';                             col = [0 0 0];     LineStyle = '-';
        otherwise,   label = field;                               col = [0 0 0];     LineStyle = '-';
    end

    plotLineWithFill(y_iv.(av_name),y_iv.(lb_name),y_iv.(ub_name), label, ...
        Color=col, FaceAlpha = 0.25, LineStyle=LineStyle, LineWidth=2,...
         x=X_all, XScale=XScale,YScale=YScale);
    hold on;

    y_entries(kIV) = struct('Color',col, 'Alpha',0.25, ...
                    'LineStyle',LineStyle, 'LineWidth',2, ...
                    'Text',label);
end
set_ylabel_2_2();
set_xlabel_2();
grid on;
axLeg2(2) = customLegend(y_entries,ax2(2),Location=LegsLoc(2),cols=LegsCols(2),FontSize=opts.FS_Legend);
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