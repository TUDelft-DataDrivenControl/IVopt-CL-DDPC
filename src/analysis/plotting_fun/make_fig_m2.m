function make_fig_m2(m2,X_all,useFillCases,noPlotCases,opts)
%MAKE_FIG_M2  Plot figure for Lf error types across cases with custom options.
%
%   make_fig_m2(m2, X_all, useFillCases, noPlotCases, opts)
%
%   Plots median and percentile bounds for different cases and Lf error types.
%   Uses custom colors, line styles, and legend placement. Optionally fills bounds.
%
%   Inputs:
%     m2           Struct with fields for each Lf error type and case, containing median and percentiles.
%     X_all        Vector of N or Re values (x-axis).
%     useFillCases Cell array of case names for which to fill percentile bounds.
%     noPlotCases  Cell array of case names to exclude from plotting.
%     opts         Options struct with fields:
%                   PlotMode      - x-axis label ('N' or 'Re') (default: 'N')
%                   fontSize      - Font size for labels (default: 15)
%                   CrameriColors - Name of Crameri color map (default: 'roma')
%                   FigPos        - Figure position [x y w h] (default: [2000 700 1000 600])
%                   fillAlpha     - Alpha for fill (default: 0.25)
%                   XScale        - X axis scale: 'log' or 'linear' (default: 'log')
%                   YScale        - Y axis scale: 'log' or 'linear' (default: 'linear')
%                   LineWidth     - Line width (default: 2)
%                   LegLoc        - Legend location (default: 'northeast')
%                   LegXY         - Legend position vector (default: [420 370])
%
%   This function requires the crameri and customLegend utilities.
arguments
    m2 struct
    X_all (:,1) double
    useFillCases cell
    noPlotCases cell
    opts.PlotMode (1,:) char {mustBeMember(opts.PlotMode,{'N','Re'})} = 'N'
    opts.fontSize (1,1) double {mustBeReal,mustBeFinite,mustBePositive} = 15
    opts.CrameriColors (1,:) char = 'roma'  % passed to crameri()
    opts.FigPos (1,4) double {mustBeReal,mustBeFinite} = [2000 700 1000 600]
    opts.fillAlpha (1,1) double {mustBeGreaterThanOrEqual(opts.fillAlpha,0),...
                                  mustBeLessThanOrEqual(opts.fillAlpha,1)} = 0.25
    opts.XScale (1,:) char {mustBeMember(opts.XScale,["log","linear"])} = 'log'
    opts.YScale (1,:) char {mustBeMember(opts.YScale,["log","linear"])} = 'linear'
    opts.LineWidth (1,1) double {mustBeReal,mustBeFinite,mustBePositive} = 2
    opts.LegLoc (1,1) string = "northeast"   % passed to customLegend
    opts.LegXY (1,2) double {mustBeReal,mustBeFinite} = [420 370] % position legend
end

% Unpack opts into variables
fontSize  = opts.fontSize;      % font size of x & y labels
cColors   = opts.CrameriColors;
Position  = opts.FigPos;
fillAlpha = opts.fillAlpha;
XScale    = opts.XScale;
YScale    = opts.YScale;
LineWidth = opts.LineWidth;
LegLoc    = opts.LegLoc;
LegXY     = opts.LegXY;
if strcmp(opts.PlotMode,'N')
    xAxisLabel = '$N$';
else
    xAxisLabel = '$\mathrm{Var}(e_k)$';
end


% --------------------- set colours and make figure -----------------------
try
    cCram = crameri(cColors, 7); % colors
catch ME
    msg = ME.message;
    msg = sprintf('Crameri color palette not recognized, resulting in:\n%s',msg);
    ME.message = msg;
   rethrow(ME)
end

fig3 = figure();
fig3.Position = Position;
tl3 = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
LfErrorTypes = fieldnames(m2);
Cases = fieldnames(m2.Up);
Cases = setdiff(Cases,noPlotCases); % only plot cases not in noPlotCases
nCases = numel(Cases);

for kEt = 1:numel(LfErrorTypes)
    Et = LfErrorTypes{kEt};
    ax3(kEt) = nexttile(tl3);

    for kC = 1:nCases
        CaseName = Cases{kC};
        
        switch CaseName
            case 'iv1',                  col = [1 0 0];    LineStyle = '-'; % red
            case 'iv2a',                 col = [0 0 0]; % black
            case 'iv2b',                 col = cCram(3,:);
            case 'iv2c',                 col = cCram(7,:);
            case {'iv3a','iv3c'},        col = cCram(1,:);
            case {'iv4a','iv4b','iv4c'}, col = cCram(2,:);
            case {'iv5a','iv5b','iv5c'}, col = cCram(4,:);
            case {'iv6a','iv6c'},        col = cCram(5,:);
            case 'CLSPC', col = cCram(6,:); LineStyle = '-'; label = 'CL-SPC';
            case 'actLf', col = [0 1 0 ];   LineStyle = '-'; label = 'actual $L_f$'; % green
        end
        if startsWith(CaseName,'iv')
            label = sprintf('$j=%s$',CaseName(3:end));
            if endsWith(CaseName,'a')
                LineStyle = '-.';
            elseif endsWith(CaseName,'b')
                LineStyle = '--';
            elseif endsWith(CaseName,'c')
                LineStyle = ':';
            end
        end
        
        % determine whether to use fill for bounds
        if any(ismember(CaseName,useFillCases))
            useFill = true;
        else
            useFill = false;
        end

        plotLineWithFill(m2.(Et).(CaseName).median,...
               m2.(Et).(CaseName).pctiles(:,6),...  25th percentile -> 6
               m2.(Et).(CaseName).pctiles(:,16),... 75th percentile -> 16
               label, Color=col, FaceAlpha = fillAlpha, LineStyle=LineStyle, LineWidth=LineWidth,...
               x=X_all, XScale=XScale,YScale=YScale,useFill = useFill);
        hold on;
        
        if kEt == 1
            entries3(kC) = struct('Color',col, 'Alpha',fillAlpha*useFill, ...
                    'LineStyle',LineStyle, 'LineWidth',LineWidth, ...
                    'Text',label);
        end

    end
    grid on;
    Et_str = append(Et(1),'_',Et(2));
    title_str = append('$i=',Et_str,'$');
    title(title_str,'Interpreter','latex','FontSize',fontSize);
    xlabel(xAxisLabel,'FontSize',fontSize,'Interpreter','latex');
    if kEt == 1
        ylabel('$\left\| \Delta L_f^{i}\right\|_2$','FontSize',fontSize,'Interpreter','latex');
    end
end
linkaxes([ax3(:)],'x');
axLeg_3 = customLegend(entries3,ax3(1),cols=2,Location=LegLoc,RelScaling=false);
axLeg_3.Position(1:2) = LegXY;
end