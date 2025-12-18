%% This script plots all of the figures that feature in the publication
% accompanying this repository
clear;
clc;
close all;

pdir = load('pdir.mat').pdir;
cd(pdir);

%% Figure 1 - example of Monte-Carlo simulation

% settings for plot_run
data_type = 'Re'; % 'N', 'Re', or 'p' represented by X below
iX = 15;          % index of data_type
ks = 5;           % seed index
noPlotCases  = {'iv2b','iv2c','iv3c','iv4b','iv4c','iv5b','iv5c','iv6c'};
subdir1 = 'Re_1e-04_1e-01_20_p_20_N_1e03_f_20_20251028_1223';
subdir1 = fullfile(pdir,'data','raw','sys1','ref0_prbs','dRe',subdir1);

plot_run; % create plot of Monte-Carlo simulation

% 
lgd = findobj(gcf, 'Type', 'Legend');
lgd.NumColumns = 3;
lgd.Location = 'southwest';
lgd.FontSize = 12;

fig1 = gcf;
fig1.Position = [680   283   940   595];

allAxes = findall(gcf, 'Type', 'axes');

axisFontSize = 14;
% Loop through and set label font sizes
for i = 1:length(allAxes)
    allAxes(i).FontSize = axisFontSize-2; % numbers on axis
    allAxes(i).XLabel.FontSize = axisFontSize;
    allAxes(i).YLabel.FontSize = axisFontSize;
end

xlim(gca,[0 length(yr0)+length(yr1)-1])

% exporting figure
set(fig1, 'Color', 'w');
exportgraphics(fig1, fullfile(pdir,'results','fig1_exampleMC.pdf'), ...
    'BackgroundColor', 'white', ...
    'ContentType', 'vector', ...
    'Resolution', 300);

%% Figure 2 - example of IVs
clearvars -except pdir
cd(pdir);

% ------------------------- settings --------------------------------------
% subdir1 = 'Re_1e-04_1e-01_20_p_20_N_1e03_f_20_20251028_1223';
% only 1 seed per Re setting, so useful for example w/ particular y0, e0, u0:
subdir1 = 'Re_1e-05_1e00_10_p_20_N_1e03_f_20_20251029_1654';
fig2_dir = fullfile(pdir,'data','raw','sys1','ref0_prbs','dRe',subdir1);
fig2_file = fullfile(fig2_dir,'processed_data.mat');
pX = 10; % index of Re in Re_all to plot this for
useFillCases = {'iv2a','iv2b'};
noPlotCases  = {'iv1'};

FS_Tick   = 8;
FS_Label  = 9;
Fs_Legend = 7;

% loading data
load(fig2_file);
% plot figure
[fig2,ax2,axLeg2] = make_fig_m0(m0, pX, useFillCases, noPlotCases,... % plots median & 25th/75th percentile over all seeds
                     LegCols=[3,3],CrameriColors='roma',... -> use subdir w/ only 1 seed per setting for example
                     LegIVnumOnly=true,FigPos=[0.05 0.05 252 300],FigUnits='points',...[-2400 -200 1000 900]
                     FS_Label=FS_Label,FS_Legend=Fs_Legend,LineWidth=1,...
                     Leg_BoxWidth=30); % 15, 12
for k=1:2
    ax2(k).FontSize = FS_Tick;
    ax2(k).YLabel.FontSize = FS_Label;
end
ax2(k).XLabel.FontSize = FS_Label;
fig2.OuterPosition(3:4) = [252 252/fig2.OuterPosition(3)*fig2.OuterPosition(4)];

% getting data from directory
dirs = dir(fig2_dir);
names = {dirs.name};
pattern = sprintf('^0*%d_',pX);
matches = ~cellfun(@isempty, regexp(names, pattern, 'once')); %match pattern
firstMatchIdx = find(matches, 1, 'first'); % Find the first match
pX_dir = fullfile(fig2_dir,names{firstMatchIdx});
cd(pX_dir);
seed_files = dir('seed_*.mat');
cd(pdir);
load(fullfile(pX_dir,seed_files(1).name),'opts','y0','e0','u0'); % makes sense, also only one seed per Re in subdir1

% overlaying e_k & y_k on plots
plot(ax2(1),e0(opts.p+1:end),'r-','DisplayName','$e_k$','LineWidth',1);
plot(ax2(1),y0(opts.p+1:end),'-', 'DisplayName','$y_k$','LineWidth',1,'Color',0.5*ones(1,3));
ax2(1).YLim = [-12.75   17];
legax2_1= legend(ax2(1),'Interpreter','latex','FontSize',Fs_Legend,...
                'IconColumnWidth',15,...
                'Position',[0.72    0.92    0.1239    0.0672]);

% overlaying u_k on plots
plot(ax2(2),u0(opts.p+1:end),'-','DisplayName','$u_k \Rightarrow U_{\mathrm{f}}$','LineWidth',1,'Color',0.5*ones(1,3));
ax2(2).YLim = [-3.25    6];
legax2_2 = legend(ax2(2),'Interpreter','latex','FontSize',Fs_Legend,...
                'IconColumnWidth',15,...
                'Position',[0.1235    0.4558    0.2057    0.0380]);

% shifting p steps forwards -> 'future' IO data
offset = opts.p; % shift right by p timesteps
offsetXAxis(ax2(1),offset)
offsetXAxis(ax2(2),offset)

xlim(ax2(1),[940 980]);

axLeg2(1).Position(1:2) = [40.3333  356.3083];
axLeg2(2).Position(1:2) = [112.3333  142.4892];

% exporting figure
set(fig2, 'Color', 'w');
exportgraphics(fig2, fullfile(pdir,'results','fig2_exampleIVs.pdf'), ...
    'BackgroundColor', 'white', ...
    'ContentType', 'vector', ...
    'Resolution', 600);

%% Figure 3
clearvars -except pdir
cd(pdir);
% ------------------------- settings --------------------------------------
data_type = 'Re';
% subdir1 = 'Re_1e-05_1e-02_10_p_20_N_1e03_f_20_20251030_1517';
subdir1 = 'Re_1e-05_1e-01_15_p_20_N_1e03_f_20_20251030_2331';
fig3_dir = fullfile(pdir,'data','raw','sys1','ref0_prbs',['d',data_type],subdir1);
fig3_file = fullfile(fig3_dir,'processed_data.mat');
load(fig3_file);
load(fullfile(fig3_dir,sprintf('d%s_settings.mat',data_type)));
[p,f,nu,ny,N] = deal(opts.p,opts.f,opts.nu,opts.ny,opts.N);

FS_Tick   = 9;
FS_Label  = 10;
Fs_Legend = 8;

[fig3,ax3,axLeg3] = make_fig_m1(m1,'Re',Re_all,N=N, YScale='log', LegCols=[2 2],...
    LegLocations=["northwest","west"],fontSize=FS_Label,FS_Legend=Fs_Legend,...
    FigPos=[50 50 252 400],Units='points');
ax3(1).YLim(2) = 2;
fig3.Units = 'points';

for k=1:2
    ax3(k).FontSize = FS_Tick;
    ax3(k).YLabel.FontSize = FS_Label*1.5;
end
ax3(k).XLabel.FontSize = FS_Label*1.5;

set(fig3, 'Color', 'w');
exportgraphics(fig3, fullfile(pdir,'results','fig3_IVdiff_dRe.pdf'), ...
    'BackgroundColor', 'white', ...
    'ContentType', 'vector', ...
    'Resolution', 600);

%% Figure 4
clearvars -except pdir
cd(pdir);
close all;
% ------------------------- settings --------------------------------------
data_type = 'Re';
% subdir1 = 'Re_1e-05_1e-02_10_p_20_N_1e03_f_20_20251030_1517';
subdir1 = 'Re_1e-05_1e-01_15_p_20_N_1e03_f_20_20251030_2331';
fig4_dir = fullfile(pdir,'data','raw','sys1','ref0_prbs',['d',data_type],subdir1);
fig4_file = fullfile(fig4_dir,'processed_data.mat');
load(fig4_file);
load(fullfile(fig4_dir,sprintf('d%s_settings.mat',data_type)));
[p,f,nu,ny,N] = deal(opts.p,opts.f,opts.nu,opts.ny,opts.N);
iX = find(Re_all >= 5e-2,1,'first');
iX_str = sprintf('iX%d',iX);
Cases = fieldnames(m4.yfhat);
ivAs = Cases(endsWith(Cases,'a') & startsWith(Cases,'iv'));
ivBs = Cases(endsWith(Cases,'b') & startsWith(Cases,'iv'));
ivCs = Cases(endsWith(Cases,'c') & startsWith(Cases,'iv'));
noPlotCases = {'actLf'};
% noPlotCases = [{'actLf'};ivBs;ivCs;{'iv1';'iv6a';'CLSPC';'TrPred'}];
% noPlotCases = [{'actLf'};ivCs; {'iv1';'CLSPC';'TrPred'}];
% noPlotCases = setdiff(noPlotCases,{'iv2a','iv2b','iv2c'});
Cases = setdiff(Cases,noPlotCases);

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
fig4 = figure('Units', 'inches', 'Position', [1, 1, 7, 9]);
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

set(fig4, 'Color', 'w');
exportgraphics(fig4, fullfile(pdir,'results','fig4_rel_yfpred_error.pdf'), ...
    'BackgroundColor', 'white', ...
    'ContentType', 'vector', ...
    'Resolution', 600);

%% Helper functions

function offsetXAxis(ax, x_offset)
%OFFSETXAXIS Shifts all X data (lines + fills) on an axis by a given offset
%
%   offsetXAxis(gca, 10)   % shifts everything right by +10 units

    % --- Shift Line objects ---
    hLines = findall(ax, 'Type', 'Line');
    for h = hLines'
        h.XData = h.XData + x_offset;
    end

    % --- Shift Patch / Fill objects ---
    hPatches = findall(ax, 'Type', 'Patch');
    for h = hPatches'
        % Case 1: has XData (fill created by 'fill' or 'area')
        if isprop(h, 'XData') && ~isempty(h.XData)
            h.XData = h.XData + x_offset;

        % Case 2: generic patch with Vertices field
        elseif isprop(h, 'Vertices') && ~isempty(h.Vertices)
            verts = h.Vertices;
            verts(:,1) = verts(:,1) + x_offset; % shift x-column
            h.Vertices = verts;
        end
    end

    % Optional: update axis limits
    ax.XLim = ax.XLim + x_offset;
end
