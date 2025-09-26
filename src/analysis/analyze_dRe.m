clear; 
close all;

%% navigate to data\raw\dRe\<subdir1>
load("pdir.mat",'pdir'); % load path of project directory
src_dir = fullfile(pdir,'src');
raw_dir = fullfile(pdir,'data','raw');
dRe_dir = fullfile(pdir,'data','raw','dRe');

% add relevant directories to path
addpath(genpath(src_dir));

% find <subdir1> candidates:
subdir1s = find_named_subdirs(dRe_dir);

% choose <subdir1> to use
subdir1 = choose_named_subdir(subdir1s);
subdir1 = fullfile(dRe_dir,subdir1); % set path to subdir1
addpath(genpath(subdir1));

%% load data for choice of dRe trials (as specified by <subdir1>)
cd(subdir1); % move to subdir1

% load dRe settings from data\raw\dRe\<subdir1>\dRe_settings.mat
load('dRe_settings.mat');
spX = spRe;
[p,f,nu,ny,N] = deal(opts.p,opts.f,opts.nu,opts.ny,opts.N);

if isfile('processed_data.mat')
    load("processed_data.mat");
else
    fprintf('Processed data file not found\n')
    fprintf('Processing data in directory\n');
    process_dRe;
end


%% Get monitor positions
monitors = get(0, 'MonitorPositions');
monPos = monitors(end,:);

%% Plotting - figure 1: example of Uf_iv (m0)

pX = 8; % index of X to plot this for
make_fig_m0(m0, pX, opts);
fig1 = gcf;
fig1.OuterPosition(1:2) = monPos(1:2) + [50 50];

%% Plotting - figure 2: quality of optimal IV approx. vs. Re (m1)
make_fig_m1(m1,N,Re_all=Re_all,YScale='linear',FigPos=fig1.Position);

%% Plotting - figure 3: ID error (m2)
useFillCases = {'iv1','CLSPC','iv2a','iv2b'};
noPlotCases = {}; 
make_fig_m2(m2, Re_all, useFillCases, noPlotCases,PlotMode='Re',FigPos=fig1.Position);

%% Plotting - figure 4: DDPC performance (m3)
useFillCases = {'iv1','CLSPC','actLf'};
noPlotCases = {}; 
make_fig_m3(m3, Re_all, useFillCases, noPlotCases,PlotMode='Re',FigPos=fig1.Position);

%% remove data path again
cd(src_dir);
rmpath(genpath(subdir1));

%% Helper functions
function subdirs = find_named_subdirs(parentDir)
    d = dir(parentDir);
    d = d([d.isdir]);
    names = {d.name};
    names = names(~ismember(names,{'.','..'}));
    pattern = ['^Re_[-0-9ep]+_[-0-9ep]+_\d+_sys_\d+_p_\d+_f_\d+_N_[-0-9ep]+_Ncl_[-0-9ep]+_Qk_[-0-9ep]+_Rk_[-0-9ep]+_dRk_[-0-9ep]+$'];
    isMatch = cellfun(@(x) ~isempty(regexp(x, pattern, 'once')), names);
    subdirs = names(isMatch);
end

function chosenDir = choose_named_subdir(subdirs)
    if numel(subdirs) > 1
        while true
            fprintf('Available subdirectories:\n');
            for i = 1:numel(subdirs)
                fprintf('  [%d] %s\n', i, subdirs{i});
            end
            choice = input('Select a subdirectory by number: ', 's');
            choiceNum = str2double(choice);
            if ~isnan(choiceNum) && choiceNum >= 1 && choiceNum <= numel(subdirs)
                break;
            else
                clc;
                fprintf('Invalid choice. Please select a number between 1 and %d.\n\n', numel(subdirs));
            end
        end
        clc;
        chosenDir = subdirs{choiceNum};
    else
        chosenDir = subdirs{1};
    end
end
