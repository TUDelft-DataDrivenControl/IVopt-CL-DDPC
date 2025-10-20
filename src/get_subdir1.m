function [subdir1,src_dir] = get_subdir1(data_type)

%% navigate to data\raw\sys#\dX\<subdir1>
load("pdir.mat",'pdir'); % load path of project directory
src_dir = fullfile(pdir,'src');
raw_dir = fullfile(pdir,'data','raw');

% choose system
sys_dirs = choose_system(raw_dir);
sys_dir  = choose_named_subdir(sys_dirs,'Choose system for which to analyze data.');
sys_dir  = fullfile(raw_dir,sys_dir);

dX_dir = fullfile(sys_dir,['d',data_type]);

% add relevant directories to path
addpath(genpath(src_dir));

% find <subdir1> candidates:
% -> subdirectories in data\raw\sys#\dX that match the naming convention from name_subdir1 in main_dX
subdir1s = find_named_subdirs_dX(data_type,dX_dir);

% choose <subdir1> to use
subdir1 = choose_named_subdir(subdir1s);
subdir1 = fullfile(dX_dir,subdir1); % set path to subdir1
addpath(genpath(subdir1));

end

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
            pattern = '^N_[-0-9ep]+_[-0-9ep]+_\d+_Re_[-0-9ep]+_p_\d+';
        case 'Re'
            pattern = '^Re_[-0-9ep]+_[-0-9ep]+_\d+_p_\d+_N_[-0-9ep]+';
        case 'p'
            pattern = '^p_\d+_\d+_\d+_Re_[-0-9ep]+_N_[-0-9ep]+';
        otherwise
            error("Data type not recognized. Choose either 'N', 'Re', or 'p'.")
    end

    
    % Keep only those that match
    isMatch = cellfun(@(x) ~isempty(regexp(x, pattern, 'once')), names);
    subdirs = names(isMatch);
end

function subdirs = choose_system(parentDir)
    % List all directories in parentDir
    d = dir(parentDir);
    d = d([d.isdir]);             % keep only directories
    names = {d.name};
    
    % Remove '.' and '..'
    names = names(~ismember(names,{'.','..'}));
    
    % Regex pattern that matches the naming convention
    pattern = '^sys\d+';
    
    % Keep only those that match
    isMatch = cellfun(@(x) ~isempty(regexp(x, pattern, 'once')), names);
    subdirs = names(isMatch);
end

function chosenDir = choose_named_subdir(subdirs,opt_str)
% CHOOSE_NAMED_SUBDIR chooses a subdirectory from subdirs
    
if numel(subdirs) > 1
    while true
        % Display options
        if nargin > 1
            fprintf([opt_str,'\n']);
        end
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
elseif ~isempty(subdirs)
    chosenDir = subdirs{1};
elseif isempty(subdirs)
    error('No subdirectories to choose from');
end

end