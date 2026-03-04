%% Choose settings for which to show simulation results
if ~exist('data_type','var')
    data_type = 'Re'; % 'N', 'Re', or 'p' represented by X below
end

if ~exist('iX','var')
    iX = 1; % index of data_type
end
if ~exist('ks','var')
    ks = 1; % seed index
end

% cases not to be plotted (will be augmented with those that are unstable)
if ~exist('noPlotCases','var')
    noPlotCases  = {'iv2b','iv2c','iv3c','iv4b','iv4c','iv5b','iv5c','iv6c'};
end

if ~exist('subdir1','var')
    subdir1 = [];
end

%% navigate to subdir1
if isempty(subdir1)
    [subdir1,~] = get_subdir1(data_type);
end
cd(subdir1);

%% load data
load(sprintf('d%s_settings.mat',data_type));

subdir2s = dir(subdir1);
isub = [subdir2s(:).isdir]; 
subdir2s = {subdir2s(isub).name};
subdir2s = subdir2s(~ismember(subdir2s,{'.','..','mfiles'}));
switch data_type
    case {'N','p'}
        if strcmp(data_type,'N')
            N = N_all(iX); X = N;
        else % -> 'p'
            p = p_all(iX); X = p;
        end
        subdir2 = choose_subdir_by_number(subdir2s, X);
        cd(subdir2); % navigate into <subdir2>

        % load file with settings in <subdir2>
        settingsFile = find_settingsFile();
        load(settingsFile);
    
    case 'Re'
        subdir2 = choose_subdir_by_iX(subdir2s, iX);
        cd(subdir2); % navigate into <subdir2>
end
load(sprintf('seed_%d.mat',seeds(ks,iX)));

f = opts.f;
Ncl = opts.Ncl;

%% Plotting

% determine unstable cases, do not plot those
fns = fieldnames(Tcl);
for k = 1:numel(fns)
    fn = fns{k};
    if ~isstable(Tcl.(fn))
        noPlotCases = [noPlotCases,{fn}];
    end
end
noPlotCases = unique(noPlotCases);

% select cases to plot
Cases  = setdiff(Cases,noPlotCases);
nCases = numel(Cases);

% determine colours to use
cCram = crameri('batlow', nCases); % colors

% time steps
Tsteps0 = 0:Nbar-1;
Tsteps1 = Nbar:Nbar+Ncl-1;
Tsteps  = [Tsteps0 Tsteps1];

% make figure and tiledlayout
figure()
tl = tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
% Define line styles and marker styles to cycle through
line_styles = {'-'};%, '--', '-.'};
marker_styles = {'o', 's', 'd', '^', 'none','v', '>', '<', 'p', 'h'};
if ~exist('marker_interval','var')
    marker_interval = 2; % Show marker every 2nd data point
end
% ------------------------------- outputs ----------------------------
ax11 = nexttile;
stairs(Tsteps,[e0 e1],'r-','DisplayName','$e_k$'); hold on; % innovation noise
stairs(Tsteps0,y0,'b-','DisplayName','past data');
for k = 1:nCases
    CaseName = Cases{k};
    switch CaseName
        case 'CLSPC'
            DispName = 'CL-SPC';
        case 'TrPred'
            DispName = 'TP';
        case 'actLf'
            DispName = 'actual $L_f$';
        otherwise
            if startsWith(CaseName,'iv')
                DispName = ['IV',CaseName(3:end)];
            else
                DispName = CaseName;
            end
    end
    
    % Cycle through line & marker styles, & marker indices
    line_style = line_styles{mod(k-1, numel(line_styles)) + 1};
    marker_style = marker_styles{mod(k-1, numel(marker_styles)) + 1};
    marker_indices = floor((k-1)/nCases*marker_interval) + 1:marker_interval:numel(Tsteps1);
    if max(marker_indices) > numel(Tsteps1)
        marker_indices = marker_indices(marker_indices <= numel(Tsteps1));
    end

    % Add markers at specified intervals
    hold on;
    plot(nan, nan, ...      just for legend
        'Color', cCram(k,:), ...
        'LineStyle', line_style, ...
        'Marker', marker_style, ...
        'MarkerSize', 6,...
        'DisplayName',DispName);
    
    stairs(Tsteps1, y_cl.(CaseName), ... plotting of stairs
        'Color', cCram(k,:), ...
        'LineStyle', line_style, ...
        'HandleVisibility', 'off');

    plot(Tsteps1(marker_indices), y_cl.(CaseName)(marker_indices), ...
        'Color', cCram(k,:), ...
        'LineStyle', 'none', ...
        'Marker', marker_style, ...
        'MarkerSize', 6,...
        'HandleVisibility', 'off');
end
stairs(0:Nbar+Ncl-1+f,[yr0 yr1],'k-','LineWidth',1.5,'DisplayName','ref.');  % references

xline(Nbar-0.5,'LineWidth',2,'Color','k','LineStyle','-.','HandleVisibility','off');
grid on;
ylabel('$y_k$','Interpreter','latex','FontSize',12);
lgd = legend;
lgd.Interpreter = 'latex';

% ------------------------------- inputs -----------------------------
ax12 = nexttile;
stairs(Tsteps0,u0,'b-'); hold on;

% YLim = ax11.YLim;
for k = 1:nCases
    CaseName = Cases{k};
    
    % Cycle through line & marker styles, & marker indices
    line_style = line_styles{mod(k-1, numel(line_styles)) + 1};
    marker_style = marker_styles{mod(k-1, numel(marker_styles)) + 1};
    marker_indices = floor((k-1)/nCases*marker_interval) + 1:marker_interval:numel(Tsteps1);
    if max(marker_indices) > numel(Tsteps1)
        marker_indices = marker_indices(marker_indices <= numel(Tsteps1));
    end

    % Add markers at specified intervals
    hold on;   
    stairs(Tsteps1, u_cl.(CaseName), ... plotting of stairs
        'Color', cCram(k,:), ...
        'LineStyle', line_style, ...
        'HandleVisibility', 'off');

    plot(Tsteps1(marker_indices), u_cl.(CaseName)(marker_indices), ...
        'Color', cCram(k,:), ...
        'LineStyle', 'none', ...
        'Marker', marker_style, ...
        'MarkerSize', 6,...
        'HandleVisibility', 'off');
end
% ylim(YLim);
stairs(Nbar:Nbar+Ncl-1+f,ur1,'k-','LineWidth',1.5);  % references

xline(Nbar-0.5,'LineWidth',2,'Color','k','LineStyle','-.');
grid on;
ylabel('$u_k$','Interpreter','latex','FontSize',12);
xlabel('Time step','Interpreter','latex','FontSize',12);

linkaxes([ax11 ax12],'x');

%% Helper functions
function chosenDir = choose_subdir_by_number(subdirs, targetNum)
% CHOOSE_SUBDIR_BY_NUMBER selects the subdirectory whose trailing integer
% matches targetNum.
%
% Inputs:
%   subdirs   - cell array of subdirectory names (strings)
%   targetNum - integer to match at the end of the subdirectory name
%
% Output:
%   chosenDir - matching subdirectory name (string). Empty if no match.

    chosenDir = '';  % default (if no match)

    for i = 1:numel(subdirs)
        % Extract trailing number using regexp
        tokens = regexp(subdirs{i}, '(\d+)$', 'tokens');
        if ~isempty(tokens)
            num = str2double(tokens{1}{1});
            if num == targetNum
                chosenDir = subdirs{i};
                return;  % stop after first match
            end
        end
    end

    if isempty(chosenDir)
        warning('No subdirectory ends with the number %d.', targetNum);
    end
end

function settingsFile = find_settingsFile(dirPath)
% FIND_MAT_FILES returns all .mat files in the specified directory ending
% with 'settings.mat'.
%
% Inputs:
%   dirPath - path to the directory (string). Defaults to pwd.
%
% Outputs:
%   settingsFiles - cell array of full paths to .mat files ending with 'settings.mat'

    if nargin < 1
        dirPath = pwd;  % default to current directory
    end

    % Get all .mat files
    settingsFile = dir(fullfile(dirPath, '*settings.mat'));

    % convert to cell array
    settingsFile = {settingsFile.name};

    if numel(settingsFile) > 1
        error('only expecting one settings file')
    end
    settingsFile = settingsFile{1};
end

function chosenDir = choose_subdir_by_iX(subdirs, iX)
% CHOOSE_SUBDIR_BY_IX selects the subdirectory whose leading zero-padded number matches iX.
% Example: for iX=3, matches '003_someName' if present.
    chosenDir = '';
    for i = 1:numel(subdirs)
        dirName = subdirs{i};
        tokens = regexp(dirName, '^(\d+)', 'tokens');
        if ~isempty(tokens)
            dirNum = str2double(tokens{1}{1});
            if dirNum == iX
                chosenDir = dirName;
                return;
            end
        end
    end
    if isempty(chosenDir)
        warning('No subdirectory found starting with number %d.', iX);
    end
end