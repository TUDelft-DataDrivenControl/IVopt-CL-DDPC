function fig = Fig_prediction_quality(fig_dir,iX,data_type,varargin)
% fig = Fig_prediction_quality(fig_dir,iX,data_type,Cases,noPlotCases)
narginchk(3,5); % at least 3, at most 5 arguments

nVararg = length(varargin);
switch nVararg
    case 1
        Cases = varargin{1};
        noPlotCases = {};
    case 2
        Cases = varargin{1};
        noPlotCases = varargin{2};
        Cases = setdiff(Cases,noPlotCases);
end

fig_file = fullfile(fig_dir,'processed_data.mat');
load(fig_file);
load(fullfile(fig_dir,sprintf('d%s_settings.mat',data_type)));
[p,f,nu,ny,N] = deal(opts.p,opts.f,opts.nu,opts.ny,opts.N);
% iX = find(Re_all >= 5e-2,1,'first');
iX_str = sprintf('iX%d',iX);
if isempty(Cases)
    Cases = fieldnames(m4.yfhat);
    ivAs = Cases(endsWith(Cases,'a') & startsWith(Cases,'iv'));
    ivBs = Cases(endsWith(Cases,'b') & startsWith(Cases,'iv'));
    ivCs = Cases(endsWith(Cases,'c') & startsWith(Cases,'iv'));
    noPlotCases = {'actLf'};
    % noPlotCases = [{'actLf'};ivBs;ivCs;{'iv1';'iv6a';'CLSPC';'TrPred'}];
    % noPlotCases = [{'actLf'};ivCs; {'iv1';'CLSPC';'TrPred'}];
    % noPlotCases = setdiff(noPlotCases,{'iv2a','iv2b','iv2c'});
    % Cases = setdiff(Cases,noPlotCases);
    % Cases = {'iv1'};
    % Cases = {'iv2a','iv2b','iv2c'};
    % Cases = {'iv3a',       'iv3c'};
    % Cases = {'iv4a','iv4b','iv4c'};
    % Cases = {'iv5a','iv5b','iv5c'};
    % Cases = {'iv6a',       'iv6c'};
    Cases = {'CLSPC','TrPred','iv1','iv2a','iv4a','iv5a','iv3a','iv6a'};
end

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
fig = figure('Units', 'inches', 'Position', [1, 1, 7, 9]);
tiledlayout(2,1,'TileSpacing','tight','Padding','tight');

% Plot mean
ax4(1) = nexttile;
hold on; box on; grid on;
h = gobjects(nCases, 1);
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

end