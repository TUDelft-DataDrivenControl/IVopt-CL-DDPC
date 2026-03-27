%% This script plots all of the figures that feature in the publication
% accompanying this repository
clear;
clc;
close all;

save_figs = false;

[src_plot_figs_dir, ~  , ~] = fileparts(which(mfilename)); % find src/util directory
cd(fullfile(src_plot_figs_dir,'..','..')); pdir = pwd; % go to project directory and save the path
addpath(genpath(fullfile(pdir,'src')));          % add src & subdirectories to path
addpath(fullfile(pdir,'bin','crameri_colours')); % add path to Crameri colour maps for plotting

%% data to be plotted:
data_type = 'Re'; % 'N', 'Re', or 'p' represented by X below
subdir1 = 'Re_1e-05_1e-01_15_p_20_N_1e03_f_20_20260302_1826';             % name of subdir1 containing data for all figures
top_dir = fullfile(pdir,'data','raw','sys6_old','ref0_prbs',['d',data_type]); % parent directory of subdir1
iX = 14; % index of data_type -> determines X (Re, N, or p) below / subdir2

%% Figure: example of Monte-Carlo simulation
ks = 5;           % seed index
marker_interval = 20;
noPlotCases  = {'iv2a','iv3a','iv5a','TrPred',...
    'iv2b','iv2c','iv3c','iv4b','iv4c','iv5b','iv5c','iv6c'};

subdir = fullfile(top_dir,subdir1);

% Define y-axis limits: ylim_y1 = {[main_axis, zoomed_axis], ylim_y2 = {[main_axis], [zoomed_axis]}
ylim_y1 = {[-14.4 15.2], [8.9 13.4]}; % outputs: general axes and inset zoom
ylim_y2 = {[-16.0 14.5], [8.5 12.6]}; % inputs: general axes and inset zoom

fig1 = Fig_sim_example(data_type,iX,ks,noPlotCases,subdir,marker_interval,ylim_y1,ylim_y2);

% conditional save of figure to results folder
save_fig(save_figs,fig1,pdir,'dinkl2.pdf'); %formerly fig1_exampleMC.pdf

%% Figure: approximation of optimal IV
clearvars -except pdir save_figs data_type subdir1 top_dir iX
cd(pdir); close all;
% ------------------------- settings --------------------------------------
fig3_dir = fullfile(top_dir,subdir1);

legPosOnAx = {[0.38, 0.05], [0.53 0.05]};
legendEntriesPerColumn = {[1 2 3],[2 2]};
[axs3, fig3] = Fig_IV_approx(data_type,fig3_dir,legPosOnAx,legendEntriesPerColumn);
axs3(1).YLim(1) = 4e-2;

% conditional save of figure to results folder
save_fig(save_figs,fig3,pdir,'dinkl3.pdf'); %formerly fig3_IVdiff_dRe.pdf

%% Figure: Lf estimates
clearvars -except pdir save_figs data_type subdir1 top_dir iX
cd(pdir); close all;

Cases = {'actLf','iv1','iv6a','iv4a','iv3c','iv3a','TrPred'};

subdir1_path = fullfile(top_dir,subdir1);

cd(subdir1_path);
load("dRe_settings.mat")
load("processed_data.mat");

fig = Fig_Lf_estimates(Cases,mLf,iX,opts);
save_fig(save_figs,fig,pdir,'dinkl4.pdf'); % formerly fig2_Lf_estimates.pdf

%% Figure: yf prediction quality
clearvars -except pdir save_figs data_type subdir1 top_dir iX
cd(pdir); close all;

% ------------------------- settings --------------------------------------
fig4_dir = fullfile(top_dir,subdir1);

Cases = {'CLSPC','TrPred','actLf','iv4a','iv1','iv3c','iv3a','iv6a'};
MainYLims  = {[-0.05 0.12],[0     3.0]};    % yLims for main axes; 'free' or [ymin ymax] vector
insetYLims = {[-0.41 0.41], 'free'};    % yLim for top & bottom-axes insets
insetPos   = {[0.11, 0.61, 0.52, 0.34], ... % [X, Y, W, H] for top-axes inset
              [0.11, 0.61, 0.52, 0.34]};  % [X, Y, W, H] for bottom-axes inset
ConnectorLocation = {'south','south'};  % connector side: 'south','north','west','east'
[fig4,max_ratio] = Fig_prediction_quality(fig4_dir,iX,data_type,Cases,[],insetYLims,insetPos,MainYLims,ConnectorLocation);
% max_ratio reports the maximum ratio of (mean / std. dev.) of the prediction error for all k and selected cases

% conditional save of figure to results folder
save_fig(save_figs,fig4,pdir,'dinkl5.pdf'); % formerly fig4_rel_yfpred_error.pdf

%% helper function

% conditional save of figure to results folder
function save_fig(save_figs,fig,pdir,fig_name)
if save_figs
    set(fig, 'Color', 'w');
    exportgraphics(fig, fullfile(pdir,'results',fig_name), ...
        'BackgroundColor', 'white', ...
        'ContentType', 'vector', ...
        'Resolution', 600);
end
end

