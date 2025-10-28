%% This script plots all of the figures that feature in the publication
% accompanying this repository
clear;
clc;
close all;

pdir = load('pdir.mat').pdir;

%% Figure 1 - example of Monte-Carlo simulation
data_type = 'Re'; % 'N', 'Re', or 'p' represented by X below
iX = 15; % index of data_type
ks = 5; % seed index
noPlotCases  = {'iv2b','iv2c','iv3c','iv4b','iv4c','iv5b','iv5c','iv6c'};
subdir1 = 'Re_1e-04_1e-01_20_p_20_N_1e03_f_20_20251028_1223';
subdir1 = fullfile(pdir,'data','raw','sys1','ref0_prbs','dRe',subdir1);

plot_run;
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

set(fig1, 'Color', 'w');
exportgraphics(fig1, fullfile(pdir,'results','fig1_exampleMC.pdf'), ...
    'BackgroundColor', 'white', ...
    'ContentType', 'vector', ...
    'Resolution', 300);
