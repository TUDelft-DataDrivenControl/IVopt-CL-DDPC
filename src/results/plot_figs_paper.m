%% This script plots all of the figures that feature in the publication
% accompanying this repository
clear;
clc;
close all;

save_figs = false;

pdir = load('pdir.mat').pdir;
cd(pdir);
addpath(genpath(pwd));

%% Figure: example of Monte-Carlo simulation
% settings for plot_run
data_type = 'Re'; % 'N', 'Re', or 'p' represented by X below
iX = 14;          % index of data_type
ks = 5;           % seed index
marker_interval = 20;
noPlotCases  = {'iv2a','iv3a','iv5a','TrPred',...
    'iv2b','iv2c','iv3c','iv4b','iv4c','iv5b','iv5c','iv6c'};

% subdir1 = 'Re_1e-04_1e-01_20_p_20_N_1e03_f_20_20251028_1223';
% subdir1 = fullfile(pdir,'data','raw','sys1','ref0_prbs','dRe',subdir1);

% subdir1 = 'Re_1e-05_1e-01_15_p_20_N_1e03_f_20_20260310_1250';
% subdir1 = fullfile(pdir,'data','raw','sys7','ref0_prbs','dRe',subdir1);

subdir1 = 'Re_1e-05_1e-01_15_p_20_N_1e03_f_20_20260302_1826';
% subdir1 = 'Re_1e-05_1e-01_15_p_20_N_1e02_f_20_20260303_1429';
% subdir1 = 'Re_1e-05_1e-01_15_p_20_N_3e02_f_20_20260303_1453';
subdir1 = fullfile(pdir,'data','raw','sys6','ref0_prbs','dRe',subdir1);

% Define y-axis limits: ylim_y1 = {[main_axis, zoomed_axis], ylim_y2 = {[main_axis], [zoomed_axis]}
ylim_y1 = {[-14.4 15.2], [8.9 13.4]}; % outputs: general axes and inset zoom
ylim_y2 = {[-16.0 14.5], [8.5 12.6]}; % inputs: general axes and inset zoom

fig1 = Fig_sim_example(data_type,iX,ks,noPlotCases,subdir1,marker_interval,ylim_y1,ylim_y2);

% conditional save of figure to results folder
save_fig(save_figs,fig1,pdir,'fig1_exampleMC.pdf');

%% Figure: example of IVs
%{
clearvars -except pdir save_figs
cd(pdir);

data_case = 4;
switch data_case
    case 1
        subdir1 = 'Re_1e-05_1e00_10_p_20_N_1e03_f_20_20251029_1654';
        sysnum = 1;
        iX = 10; % index of Re in Re_all to plot this for
        xyLims = {[940 980],...     xLim
                  [-12.75   17],... yLim of top axes
                  [-3.25    6]};  % yLim of bottom axes
        LegsPos = {[0.72    0.92    0.1239    0.0672],... % pos. of smaller legend top axes
                  [40.3333  356.3083],...                 % pos. of larger  legend top axes
                  [0.1235    0.4558    0.2057    0.0380],...pos. of smaller legend bottom axes
                  [112.3333  142.4892]};                  % pos. of larger  legend bottom axes
    case 2
        subdir1 = 'Re_1e-05_1e-01_15_p_20_N_1e03_f_20_20260302_1826';
        sysnum = 6;
        iX = 14;
        xyLims = {[936 976],...
                  [-2.3  3.5],...
                  [-2    4]};
        LegsPos = {[0.72    0.92    0.1239    0.0672],... % pos. of smaller legend top axes
                  [40.3333  356.3083],...                 % pos. of larger  legend top axes
                  [0.1235    0.4558    0.2057    0.0380],...pos. of smaller legend bottom axes
                  [112.3333  142.4892]};                  % pos. of larger  legend bottom axes
    case 3
        subdir1 = 'Re_1e-05_1e-01_15_p_20_N_1e02_f_20_20260303_1429';
        sysnum = 6;
        iX = 14;
        xyLims = {[70 96.5],...     xLim
                  [-2   3],... yLim of top axes
                  [-2  3.65]};  % yLim of bottom axes
        LegsPos = {[0.72    0.92    0.1239    0.0672],... % pos. of smaller legend top axes
                  [40.3333  356.3083],...                 % pos. of larger  legend top axes
                  [0.1235    0.4558    0.2057    0.0380],...pos. of smaller legend bottom axes
                  [112.3333  142.4892]};                  % pos. of larger  legend bottom axes
    case 4
        subdir1 = 'Re_1e-05_1e-01_15_p_20_N_3e02_f_20_20260303_1453';
        sysnum = 6;
        iX = 14;
        xyLims = {[70 96.5],...     xLim
                  [-2   3],... yLim of top axes
                  [-2  3.65]};  % yLim of bottom axes
        LegsPos = {[0.72    0.92    0.1239    0.0672],... % pos. of smaller legend top axes
                  [40.3333  356.3083],...                 % pos. of larger  legend top axes
                  [0.1235    0.4558    0.2057    0.0380],...pos. of smaller legend bottom axes
                  [112.3333  142.4892]};                  % pos. of larger  legend bottom axes
end
useFillCases = {'iv2a','iv2b'};
noPlotCases  = {'iv1'};

data_type = regexp(subdir1, '^([^_]+)_', 'tokens', 'once'); data_type = data_type{1};
fig2_dir = fullfile(pdir,'data','raw',sprintf('sys%d',sysnum),'ref0_prbs',['d',data_type],subdir1);
fig2 = Fig_IV_example(fig2_dir, iX, useFillCases, noPlotCases, pdir, xyLims, LegsPos);

% conditional save of figure to results folder
save_fig(save_figs,fig2,pdir,'fig2_exampleIVs.pdf');
%}


%% Figure: approximation of optimal IV
clearvars -except pdir save_figs
cd(pdir); close all;
% ------------------------- settings --------------------------------------
data_type = 'Re';
% subdir1 = 'Re_1e-05_1e-01_15_p_20_N_1e03_f_20_20260310_1250';
% fig3_dir = fullfile(pdir,'data','raw','sys7','ref0_prbs',['d',data_type],subdir1);
% subdir1 = 'Re_1e-05_1e-01_15_p_50_N_1e03_f_20_20260310_1423';
% fig3_dir = fullfile(pdir,'data','raw','sys8','ref0_prbs',['d',data_type],subdir1);
subdir1 = 'Re_1e-05_1e-01_15_p_20_N_1e03_f_20_20260302_1826';
fig3_dir = fullfile(pdir,'data','raw','sys6','ref0_prbs',['d',data_type],subdir1);

legPosOnAx = {[0.38, 0.05], [0.53 0.05]};
legendEntriesPerColumn = {[1 2 3],[2 2]};
[axs3, fig3] = Fig_IV_approx(data_type,fig3_dir,legPosOnAx,legendEntriesPerColumn);
axs3(1).YLim(1) = 4e-2;

% conditional save of figure to results folder
save_fig(save_figs,fig3,pdir,'fig3_IVdiff_dRe.pdf');

%% Figure: Lf estimates
clearvars -except pdir save_figs
cd(pdir); close all;

% Cases = {'actLf','iv1','iv6a','iv5a','iv4a','iv3a','TrPred'};
Cases = {'actLf','iv1','iv6a','iv4a','iv3c','iv3a','TrPred'};
% subdir1 = 'Re_1e-05_1e-01_15_p_20_N_1e03_f_20_20260310_1250';
% subdir1 = fullfile(pdir,'data','raw','sys7','ref0_prbs','dRe',subdir1);
% subdir1 = 'Re_1e-05_1e-01_15_p_50_N_1e03_f_20_20260310_1423';
% subdir1 = fullfile(pdir,'data','raw','sys8','ref0_prbs','dRe',subdir1);
subdir1 = 'Re_1e-05_1e-01_15_p_20_N_1e03_f_20_20260302_1826';
subdir1 = fullfile(pdir,'data','raw','sys6','ref0_prbs','dRe',subdir1);
cd(subdir1);
load("dRe_settings.mat")
load("processed_data.mat");
iX = 14;
fig = Fig_Lf_estimates(Cases,mLf,iX,opts);
save_fig(save_figs,fig,pdir,'fig2_Lf_estimates.pdf');

%% Figure: yf prediction quality
clearvars -except pdir save_figs
cd(pdir); close all;

% ------------------------- settings --------------------------------------
data_type = 'Re';
% subdir1 = 'Re_1e-05_1e-01_15_p_20_N_1e03_f_20_20260310_1250';
% fig4_dir = fullfile(pdir,'data','raw','sys7','ref0_prbs',['d',data_type],subdir1);
% subdir1 = 'Re_1e-05_1e-01_15_p_50_N_1e03_f_20_20260310_1423';
% fig4_dir = fullfile(pdir,'data','raw','sys8','ref0_prbs',['d',data_type],subdir1);
subdir1 = 'Re_1e-05_1e-01_15_p_20_N_1e03_f_20_20260302_1826'; % <-- removed normalization in computation in m4. clearer figure!
fig4_dir = fullfile(pdir,'data','raw','sys6','ref0_prbs',['d',data_type],subdir1);
iX = 14;

Cases = {'actLf','CLSPC','TrPred','iv1','iv4a','iv3c','iv3a','iv6a'};
MainYLims  = {[-0.05 0.12],[0     3.0]};    % yLims for main axes; 'free' or [ymin ymax] vector
insetYLims = {'free', 'free'};    % yLim for top & bottom-axes insets
insetPos   = {[0.08, 0.58, 0.55, 0.37], ... % [X, Y, W, H] for top-axes inset
              [0.08, 0.58, 0.55, 0.37]};  % [X, Y, W, H] for bottom-axes inset
ConnectorLocation = {'south','south'};  % connector side: 'south','north','west','east'
[fig4,max_ratio] = Fig_prediction_quality(fig4_dir,iX,data_type,Cases,[],insetYLims,insetPos,MainYLims,ConnectorLocation);
max_ratio

% conditional save of figure to results folder
save_fig(save_figs,fig4,pdir,'fig4_rel_yfpred_error.pdf');

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

