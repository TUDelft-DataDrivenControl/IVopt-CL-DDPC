clear; 
close all;
data_type = 'p'; % 'N', 'Re', or 'p' represented by X below

%% navigate to data\raw\d<Re \ N \ p>\<subdir1>
load("pdir.mat",'pdir'); % load path of project directory
src_dir = fullfile(pdir,'src');
raw_dir = fullfile(pdir,'data','raw');
dX_dir = fullfile(pdir,'data','raw',['d',data_type]);

% add relevant directories to path
addpath(genpath(src_dir));

% find <subdir1> candidates:
% -> subdirectories in data\raw\dX that match the naming convention from name_subdir1 in main_dX
subdir1s = find_named_subdirs_dX(data_type,dX_dir);

% choose <subdir1> to use
subdir1 = choose_named_subdir(subdir1s);
subdir1 = fullfile(dX_dir,subdir1); % set path to subdir1
addpath(genpath(subdir1));

%% load data for choice of dRe \ dN trials (as specified by <subdir1>)
cd(subdir1); % move to subdir1

% load dN settings from data\raw\dN\<subdir1>\dN_settings.mat
load(sprintf('d%s_settings.mat',data_type));

switch data_type
    case 'N'
        spX = spN;
        [p,f,nu,ny] = deal(opts.p,opts.f,opts.nu,opts.ny);
        X_all = N_all;
    case 'Re'
        spX = spRe;
        [p,f,nu,ny,N] = deal(opts.p,opts.f,opts.nu,opts.ny,opts.N);
        X_all = Re_all;
    case 'p'
        spX = spP;
        [f,nu,ny,N] = deal(opts.f,opts.nu,opts.ny,opts.N);
        X_all = p_all;
end

if isfile('processed_data.mat')
    load("processed_data.mat");
else
    fprintf('Processed data file not found\n')
    fprintf('Processing data in directory\n');
    [m0,m1,m2,m3] = process_dX(data_type,seeds,X_all,opts);
end

%% Get monitor positions
monitors = get(0, 'MonitorPositions');
monPos = monitors(1,:);

%% Plotting - figure 1: example of Uf_iv (m0)
% possible cases: see ivs in CaseDefinitions.m

% ---------------------- user-defined plotting parameters -----------------
% useFillCases: cases for which to show bounds
% noPlotCases:  cases to exclude from final plot
switch data_type
    case 'N'
        pX = 5; % index of N in N_all to plot this for
        useFillCases = {'iv1','iv2a','iv2b'};
        noPlotCases  = {};

        % plot figure
        [fig1,ax1,axLeg1] = make_fig_m0(m0, pX, useFillCases, noPlotCases, LegCols=[2,2]);

    case 'Re'
        pX = 8; % index of Re in Re_all to plot this for
        useFillCases = {'iv1','iv2a','iv2b'};
        noPlotCases  = {};

        % plot figure
        [fig1,ax1,axLeg1] = make_fig_m0(m0, pX, useFillCases, noPlotCases, LegCols=[2,2]);

    case 'p'
        pX = 8; % index of p in p_all to plot this for
        useFillCases = {'iv1','iv2a','iv2b'};
        noPlotCases  = {};

        % plot figure
        [fig1,ax1,axLeg1] = make_fig_m0(m0, pX, useFillCases, noPlotCases, LegCols=[2,2]);

end
fig1.OuterPosition(1:2) = monPos(1:2) + [50 50];

%% Plotting - figure 2: quality of optimal IV approx. vs. X (m1)
% -> difference of Uf_iv w.r.t. Uf_iv2a
% -> difference of Yf_iv w.r.t. Yf_iv2b
switch data_type
    case 'N'
        [fig2,ax2,axLeg2] = make_fig_m1(m1,'N', N_all,                                     LegCols=[3 3],LegLocations="west");
    case 'Re'
        [fig2,ax2,axLeg2] = make_fig_m1(m1,'Re',Re_all,N=N,                YScale='linear',LegCols=[3 3],LegLocations=["west","northwest"]);
        axLeg2(1).Position(1:2) = axLeg2(1).Position(1:2) + [-10 -10];
        axLeg2(2).Position(1:2) = axLeg2(2).Position(1:2) + [-10 -20];
    case 'p'
        [fig2,ax2,axLeg2] = make_fig_m1(m1,'p', p_all, N=N,XScale='linear',YScale='linear',LegCols=[3 3],LegLocations=["west","northwest"]);
        axLeg2(1).Position(1:2) = axLeg2(1).Position(1:2) + [-10 -10];
        axLeg2(2).Position(1:2) = axLeg2(2).Position(1:2) + [-10 -20];
end
fig2.OuterPosition(1:2) = monPos(1:2) + [50 100];

%% Plotting - figure 3: ID error (m2)
% possible cases: see CaseDefinitions.m

% ---------------------- user-defined plotting parameters -----------------
% useFillCases: cases for which to show bounds
% noPlotCases:  cases to exclude from final plot
switch data_type
    case 'N'
        useFillCases = {'iv1','CLSPC','iv2a','iv2b'};
        noPlotCases  = {};

        % plotting figure
        [fig3,ax3,axLeg3] = make_fig_m2(m2, X_all, useFillCases, noPlotCases,PlotMode=data_type,LegCols=2);

    case 'Re'
        useFillCases = {'iv1','CLSPC','iv2a','iv2b'};
        noPlotCases  = {};

        % plotting figure
        [fig3,ax3,axLeg3] = make_fig_m2(m2, X_all, useFillCases, noPlotCases,PlotMode=data_type,LegCols=1);
        axLeg3.Position(1) = axLeg3.Position(1) - 80;

    case 'p'
        useFillCases = {'iv1','CLSPC','iv2a','iv2b'};
        noPlotCases  = {};

        % plotting figure
        [fig3,ax3,axLeg3] = make_fig_m2(m2, X_all, useFillCases, noPlotCases,XScale='linear',PlotMode=data_type,LegCols=2);
end
fig3.OuterPosition(1:2) = monPos(1:2) + [50 150]; % repositioning figure

%% Plotting - figure 4: DDPC performance (m3)
% possible cases: see CaseDefinitions.m

% ---------------------- user-defined plotting parameters -----------------
% useFillCases: cases for which to show bounds
% noPlotCases:  cases to exclude from final plot
switch data_type
    case 'N'
        useFillCases = {'iv1','CLSPC','actLf'};
        noPlotCases  = {};
        XScale = 'log';
        YScale = 'linear';
        LegCols = 2;
    case 'Re'
        useFillCases = {'iv1','CLSPC','actLf'};
        noPlotCases  = {};
        XScale = 'log';
        YScale = 'log';
        LegCols = 2;
    case 'p'
        useFillCases = {'iv1','CLSPC','actLf'};
        noPlotCases  = {'iv1'};
        XScale = 'linear';
        YScale = 'linear';
        LegCols = 2;
end

% --------------------------- plotting of figure --------------------------
[fig4,ax4,axLeg4] = make_fig_m3(m3, X_all, useFillCases, noPlotCases,PlotMode=data_type,XScale=XScale,YScale=YScale,LegCols=LegCols);
fig4.OuterPosition(1:2) = monPos(1:2) + [50 200];

%% remove data path again
cd(src_dir);
rmpath(genpath(subdir1));

%% Helper functions
function subdirs = find_named_subdirs_dX(data_type,parentDir)
% FIND_NAMED_SUBDIRS returns subdirectories in parentDir that match
% the naming convention from name_subdir1 in main_dN
%
% Example:
%   subdirs = find_named_subdirs_dX('N',pwd);

    % List all directories in parentDir
    d = dir(parentDir);
    d = d([d.isdir]);             % keep only directories
    names = {d.name};
    
    % Remove '.' and '..'
    names = names(~ismember(names,{'.','..'}));
    
    % Regex pattern that matches the naming convention
    switch data_type
        case 'N'        
            %   Example: N_1p0e0_1p0e1_50_sys_1_Re_1p0e3_p_2_f_3_Ncl_1p0e2_Qk_1p0e1_Rk_5p0e0_dRk_1p0e0
            pattern = ['^N_[-0-9ep]+_[-0-9ep]+_\d+_sys_\d+_' ...   % Nmin, Nmax, nN, sys
                    'Re_[-0-9ep]+_p_\d+_f_\d+_' ...                % Re, p, f
                    'Ncl_[-0-9ep]+_Qk_[-0-9ep]+_Rk_[-0-9ep]+_dRk_[-0-9ep]+$']; % Ncl, Qk, Rk, dRk
        case 'Re'
            pattern = '^Re_[-0-9ep]+_[-0-9ep]+_\d+_sys_\d+_p_\d+_f_\d+_N_[-0-9ep]+_Ncl_[-0-9ep]+_Qk_[-0-9ep]+_Rk_[-0-9ep]+_dRk_[-0-9ep]+$';
        case 'p'
            pattern = ['^p_\d+_\d+_\d+_sys_\d+_' ...               % pmin, pmax, nP, sys
                    'Re_[-0-9ep]+_N_[-0-9ep]+_f_\d+_' ...           % Re, N, f
                    'Ncl_[-0-9ep]+_Qk_[-0-9ep]+_Rk_[-0-9ep]+_dRk_[-0-9ep]+$']; % Ncl, Qk, Rk, dRk
        otherwise
            error("Data type not recognized. Choose either 'N', 'Re', or 'p'.")
    end

    
    % Keep only those that match
    isMatch = cellfun(@(x) ~isempty(regexp(x, pattern, 'once')), names);
    subdirs = names(isMatch);
end

function chosenDir = choose_named_subdir(subdirs)
% CHOOSE_NAMED_SUBDIR chooses a subdirectory from subdirs
    
if numel(subdirs) > 1
    while true
        % Display options
        fprintf('Available subdirectories:\n');
        for i = 1:numel(subdirs)
            fprintf('  [%d] %s\n', i, subdirs{i});
        end

        % Prompt user
        choiceNum = input('Select a subdirectory by number: ');

        % Validate choice
        if ~isnan(choiceNum) && choiceNum >= 1 && choiceNum <= numel(subdirs)
            % Valid choice -> exit loop
            break;
        else
            clc;
            fprintf('Invalid choice. Please select a number between 1 and %d.\n\n', numel(subdirs));
        end
    end
    clc;
    % Get chosen directory (full path)
    chosenDir = subdirs{choiceNum};
else
    chosenDir = subdirs{1};
end

end
