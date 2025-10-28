clear; 
close all;
data_type = 'Re'; % 'N', 'Re', or 'p' represented by X below
overwrite = false;

if ismember('SlurmProfile1',parallel.clusterProfiles) % running on cluster?
    onCluster = true;
else
    onCluster = false;
end

%% navigate to data\raw\sys#\ref0_<>\dX\<subdir1>
[subdir1,src_dir] = get_subdir1(data_type);
cd(subdir1); % move to subdir1
fprintf('Analyzing data in %s\n',subdir1(length(src_dir)-3:end));

%% load data for choice of dX trials (as specified by <subdir1>)

% load dX settings from data\raw\sys#\ref0_<>\dX\<subdir1>\dX_settings.mat
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
nX = numel(X_all);

% get data structures - load or by processing data
if isfile('processed_data.mat')
    expectedVars = {'m0','m1','m2','mLf','m3','m4'};
    availVars = {whos('-file','processed_data.mat').name};
    missingVars = setdiff(expectedVars, availVars);

    if isempty(missingVars)
        % all processed data available, load it
        fprintf('All necessary processed data found, loading file\n')
        load("processed_data.mat");

    elseif ~isempty(missingVars) && overwrite
        % overwrite existing processed data
        fprintf('Not all necessary processed data found, overwriting file\n')
        [m0,m1,m2,mLf,m3,m4] = process_dX(data_type,seeds,X_all,opts,plant);
        
    else
        % add only missing processed data
        fprintf('Not all necessary processed data found, adding missing data to file\n')
        [m0,m1,m2,mLf,m3,m4] = process_dX(data_type,seeds,X_all,opts,plant,missingVars); % available variables will be empty
        clear(availVars{:}); % clear empty available variables before loading them
        load("processed_data.mat",availVars{:}); % load existing variables
    end

else
    fprintf('Processed data file not found\n')
    fprintf('Processing data in directory\n');
    [m0,m1,m2,mLf,m3,m4] = process_dX(data_type,seeds,X_all,opts,plant);
end

if onCluster
    % do not attempt plotting if running on cluster
    return;
end

%% Preliminaries before plotting
monitors = get(0, 'MonitorPositions');
monPos = monitors(2,:);

Cases = fieldnames(m4.yfhat);
ivBs = Cases(endsWith(Cases,'b') & startsWith(Cases,'iv'));
ivCs = Cases(endsWith(Cases,'c') & startsWith(Cases,'iv'));

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
        pX = 10; % index of Re in Re_all to plot this for
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
        useFillCases = {'iv1','CLSPC','TrPred','iv2a','iv2b'};
        noPlotCases  = [ivBs;ivCs;{'actLf';'iv1'}];

        % plotting figure
        [fig3,ax3,axLeg3] = make_fig_m2(m2, X_all, p*ones(nX,1), f, useFillCases, noPlotCases,PlotMode=data_type,LegCols=2);

    case 'Re'
        useFillCases = {'iv1','CLSPC','iv2a','iv2b'};
        noPlotCases  = {'actLf'};

        % plotting figure
        [fig3,ax3,axLeg3] = make_fig_m2(m2, X_all, p*ones(nX,1), f, useFillCases, noPlotCases,PlotMode=data_type,LegCols=1);
        axLeg3.Position(1) = axLeg3.Position(1) - 80;

    case 'p'
        useFillCases = {'iv1','CLSPC','iv2a','iv2b'};
        noPlotCases  = {'actLf'};

        % plotting figure
        [fig3,ax3,axLeg3] = make_fig_m2(m2, X_all, p_all, f, useFillCases, noPlotCases,XScale='linear',PlotMode=data_type,LegCols=2);
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
        noPlotCases  = [ivBs;ivCs;{'iv1'}];
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
        useFillCases = {'CLSPC','actLf','iv2a'};
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