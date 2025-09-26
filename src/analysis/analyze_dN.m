clear; 
close all;

%% navigate to data\raw\dN\<subdir1>
load("pdir.mat",'pdir'); % load path of project directory
src_dir = fullfile(pdir,'src');
raw_dir = fullfile(pdir,'data','raw');
dN_dir = fullfile(pdir,'data','raw','dN');

% add relevant directories to path
addpath(genpath(src_dir));

% find <subdir1> candidates:
% -> subdirectories in data\raw\dN that match the naming convention from name_subdir1 in main_dN
subdir1s = find_named_subdirs(dN_dir);

% choose <subdir1> to use
subdir1 = choose_named_subdir(subdir1s);
subdir1 = fullfile(dN_dir,subdir1); % set path to subdir1
addpath(genpath(subdir1));

%% load data for choice of dN trials (as specified by <subdir1>)
cd(subdir1); % move to subdir1

% load dN settings from data\raw\dN\<subdir1>\dN_settings.mat
load('dN_settings.mat');
spX = spN;
[p,f,nu,ny] = deal(opts.p,opts.f,opts.nu,opts.ny);

if isfile('processed_data.mat')
    load("processed_data.mat");
else
    fprintf('Processed data file not found\n')
    fprintf('Processing data in directory\n');
    process_dN;
end

%% Plotting - figure 1: example of Uf_iv (m0)

pX = 5; % index of X to plot this for
opts.N = N_all(pX);
make_fig_m0(m0, pX, opts);

%% Plotting - figure 2: quality of optimal IV approx. vs. N (m1)
% -> difference of Uf_iv w.r.t. Uf_iv2a
% -> difference of Yf_iv w.r.t. Yf_iv2b
make_fig_m1(m1,N_all);

%% Plotting - figure 3: ID error (m2)
% possible cases:
%    Name   |       Description
% ----------|----------------------------------------------------------
%   iv1     | SPC using open-loop IV
%   iv2a    | SPC using optimal IV
%   iv2b    | SPC using optimal IV + Yf_iv
%   iv2c    | SPC using optimal IV + Yf_iv + 2SLS
%   iv3a    | SPC using LCF-IV
%   iv3c    | SPC using LCF-IV + 2SLS
%   iv4a    | SPC using approx. opt. IV w/o controller info.
%   iv4b    | SPC using approx. opt. IV w/o controller info. + Yf_iv
%   iv4c    | SPC using approx. opt. IV w/o controller info. + Yf_iv + 2SLS
%   iv5a    | SPC using approx. opt. IV w/  controller info.
%   iv5b    | SPC using approx. opt. IV w/  controller info. + Yf_iv
%   iv5c    | SPC using approx. opt. IV w/  controller info. + Yf_iv + 2SLS
%   iv6a    | SPC using basic IV: future reference
%   iv6c    | SPC using basic IV: future reference + 2SLS
%   CLSPC   | CL-SPC
%   actLf   | SPC using the actual matrix Lf

% ---------------------- user-defined plotting parameters -----------------

% set cases for which to show bounds
useFillCases = {'iv1','CLSPC','iv2a','iv2b'};

% set cases to exclude from final plot
noPlotCases = {}; 

% --------------------------- plotting of figure --------------------------
make_fig_m2(m2, N_all, useFillCases, noPlotCases);

%% Plotting - figure 4: DDPC performance (m3)
% possible cases:
%    Name   |       Description
% ----------|----------------------------------------------------------
%   iv1     | SPC using open-loop IV
%   iv2a    | SPC using optimal IV
%   iv2b    | SPC using optimal IV + Yf_iv
%   iv2c    | SPC using optimal IV + Yf_iv + 2SLS
%   iv3a    | SPC using LCF-IV
%   iv3c    | SPC using LCF-IV + 2SLS
%   iv4a    | SPC using approx. opt. IV w/o controller info.
%   iv4b    | SPC using approx. opt. IV w/o controller info. + Yf_iv
%   iv4c    | SPC using approx. opt. IV w/o controller info. + Yf_iv + 2SLS
%   iv5a    | SPC using approx. opt. IV w/  controller info.
%   iv5b    | SPC using approx. opt. IV w/  controller info. + Yf_iv
%   iv5c    | SPC using approx. opt. IV w/  controller info. + Yf_iv + 2SLS
%   iv6a    | SPC using basic IV: future reference
%   iv6c    | SPC using basic IV: future reference + 2SLS
%   CLSPC   | CL-SPC
%   actLf   | SPC using the actual matrix Lf

% ---------------------- user-defined plotting parameters -----------------

% set cases for which to show bounds
useFillCases = {'iv1','CLSPC','actLf'};

% set cases to exclude from final plot
noPlotCases = {}; 

% --------------------------- plotting of figure --------------------------
make_fig_m3(m3, N_all, useFillCases, noPlotCases);

%% remove data path again
cd(src_dir);
rmpath(genpath(subdir1));

%% Helper functions
function subdirs = find_named_subdirs(parentDir)
% FIND_NAMED_SUBDIRS returns subdirectories in parentDir that match
% the naming convention from name_subdir1 in main_dN
%
% Example:
%   subdirs = find_named_subdirs(pwd);

    % List all directories in parentDir
    d = dir(parentDir);
    d = d([d.isdir]);             % keep only directories
    names = {d.name};
    
    % Remove '.' and '..'
    names = names(~ismember(names,{'.','..'}));
    
    % Regex pattern that matches the naming convention
    %   Example: N_1p0e0_1p0e1_50_sys_1_Re_1p0e3_p_2_f_3_Ncl_1p0e2_Qk_1p0e1_Rk_5p0e0_dRk_1p0e0
    pattern = ['^N_[-0-9ep]+_[-0-9ep]+_\d+_sys_\d+_' ...   % Nmin, Nmax, nN, sys
               'Re_[-0-9ep]+_p_\d+_f_\d+_' ...            % Re, p, f
               'Ncl_[-0-9ep]+_Qk_[-0-9ep]+_Rk_[-0-9ep]+_dRk_[-0-9ep]+$']; % Ncl, Qk, Rk, dRk
    
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
        choice = input('Select a subdirectory by number: ', 's');

        % Validate choice
        choiceNum = str2double(choice);
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
