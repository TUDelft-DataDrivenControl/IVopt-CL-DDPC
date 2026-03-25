function [m1,mLf,m3,m4] = process_dX(data_type,seeds,X_all,opts,plant,OutVars)
arguments
    data_type (1,:) char {mustBeMember(data_type,{'N','Re','p'})}
    seeds double {mustBePositive,mustBeInteger,mustBeMatrix}
    X_all (1,:) double {mustBePositive}
    opts (1,1) struct
    plant ss
    OutVars (1,:) string {mustBeMember(OutVars,{'m1','mLf','m3','m4'})} = {'m1','mLf','m3','m4'};
end
spX = size(seeds,1);
[f,nu,ny] = deal(opts.f,opts.nu,opts.ny);
switch data_type
    case 'N'
        p   = opts.p;
    case 'Re'
        [N,p] = deal(opts.N,opts.p);
    case 'p'
        N = opts.N;
end
[A,BK,C,~] = ssdata(plant);
K = BK(:,nu+1:end);
Hf = make_blk_tril_toeplitz(A,K,C,eye(ny),f);

%% iterate over Re \ N values - initial data processing
subdir1 = pwd; % should be data/sys#/raw/dX/<subdir1>

% get all <subdir2> directories in data/raw/sys#/dX/<subdir1>
subdir2s = dir(subdir1);
isub = [subdir2s(:).isdir]; 
subdir2s = {subdir2s(isub).name};
subdir2s = subdir2s(~ismember(subdir2s,{'.','..'}));
subdir2s = subdir2s(~cellfun(@isempty, regexp(subdir2s, '^[0-9]+_')));

% iterate over all Re \ N \ p values
nX = numel(X_all);

% find Cases used
%  -> expect selection of {'iv1','iv2a','iv2b','iv2c','iv3a','iv3c','iv4a','iv4b','iv4c',...
%                          'iv5a','iv5b','iv5c','iv6a','iv6c','CLSPC','actLf','TrPred'};
seed_mat_files = dir(fullfile(pwd,subdir2s{1},'seed_*.mat'));
Cases = load(fullfile(pwd,subdir2s{1},seed_mat_files(1).name),'Cases').Cases;
num_Cases = numel(Cases);

%% initializing measures
pctiles = 0:5:100;

% --- IV definitions
Uf_ivs = {'iv1','iv2a','iv2c','iv3c','iv4a','iv4c','iv5a','iv5c','iv6c'}; Uf_ivs = intersect(Uf_ivs,Cases);
Yf_ivs = {'iv2b','iv3a','iv4b','iv5b','iv6a'};                            Yf_ivs = intersect(Yf_ivs,Cases);
num_Uf_ivs = numel(Uf_ivs); % needed for nested for loop inside parfor
num_Yf_ivs = numel(Yf_ivs); % needed for nested for loop inside parfor

% ======================== initialize measure 1 (m1) ======================
% -> how well IV approximates optimal IV
% -> structure: m1.<Uf/Yf>.<IVname>.data     (nX, spX)
%                                  .mean     (nX, 1)
%                                  .median   (nX, 1)
%                                  .pctiles  (nX, num_pctiles)
%     example: m1.Uf.iv1.data(iX,ks)
% -> sliceable containers:
%    m1_Uf_data(num_Uf_ivs, nX, spX)
%    m1_Yf_data(num_Yf_ivs, nX, spX)

% --- Preallocate "sliced" containers ---
% Combine Uf and Yf into a single array for easier handling: (num_Uf_ivs+num_Yf_ivs) x nX x spX
m1_UYf = zeros(num_Uf_ivs + num_Yf_ivs, nX, spX);
mLf_data = cell(nX,1); % to save Lf matrices. each cell of size num_Cases,spX,size(Lf,1),size(Lf,2)

% ======================== initialize measure 3 (m3) ======================
% -> DDPC performance
% -> structure: m3.<costType>.<caseName>.data    (nX, spX)
%                                       .mean    (nX, 1)
%                                       .median  (nX, 1)
%                                       .pctiles (nX, num_pctiles)
%      example: m3.cost_u.iv1.data(iX,ks)
% -> sliceable containers:
%    m3_data(num_cost_types, num_Cases, nX, spX)

% --- cost types and cases
cost_types = {'cost_u','cost_y','cost_tot'};

% --- Preallocate sliceable containers ---
m3_data = zeros(numel(cost_types), numel(Cases), nX, spX);

% ======================== measure 4 (prediction error) ===================
% -> prediction error
m4_data = zeros(nX,ny*f,num_Cases,3);

%% -------------------- loop over Re \ N \ p values -----------------------
% Check and start parallel pool if needed
if ismember('SlurmProfile1',parallel.clusterProfiles) % on cluster?
    if ~isempty(gcp('nocreate'))
        delete(gcp('nocreate'));
    end
    myCluster = parcluster('SlurmProfile1');
    nworkers = 20;
    myCluster.SubmitArguments = SlurmSubmitArgs(['calc_',data_type],15,ntasks=nworkers,cpt=1,GB=3.8);
    parpool(myCluster,nworkers);
else
    cluster = parcluster('local');
    max_workers = cluster.NumWorkers;
    if ~isempty(gcp('nocreate'))
        curr_pool = gcp('nocreate');
        if curr_pool.NumWorkers < max_workers
            delete(curr_pool);
            parpool('local',max_workers);
        end
    else
        parpool('local', max_workers);
    end
end

% start iterating over Re \ N \ p values
for iX = 1:nX
    switch data_type
        case {'N','p'}
            if strcmp(data_type,'N')
                N = X_all(iX); X = N;
            else % -> 'p'
                p = X_all(iX); X = p;
            end
            % choose subdir2 corresponding with N value
            subdir2 = choose_subdir_by_number(subdir2s, X);
            cd(subdir2); % navigate into <subdir2>

            % load file with settings in <subdir2>
            settingsFile = find_settingsFile();
            load(settingsFile);

        case 'Re'
            X = X_all(iX);
            % choose subdir2 corresponding with Re value
            subdir2 = choose_subdir_by_iX(subdir2s, iX);
            cd(subdir2); % navigate into <subdir2>
    end
    
    % ======================== processing m1, m3 ==================
    fprintf('Processing %s index [%d/%d] (%s = %g) in subdir: %s\n', data_type, iX, nX, data_type, X, subdir2);
    
    % iterates over noise realizations
    [m1_UYf(:, iX, :), mLf_data{iX},m3_data(:, :, iX, :), m4_data(iX,:,:,:,:)] ...
     = process_m_all(Uf_ivs,Yf_ivs,Cases,cost_types,iX,nu,ny,p,f,seeds,spX,Hf,OutVars);

    cd(subdir1);
end
if ismember('SlurmProfile1',parallel.clusterProfiles)
    delete(gcp('nocreate')); % close parallel pool
end

%% ============= Data processing and transforming to structures ============
[m1,mLf,m3,m4] = deal([]); % overwritten if in OutVars
for k = 1:numel(OutVars)
    OutVar = OutVars{k};
    switch OutVar
        
        % ----- m1 (quality of the approximation of the optimal IV) ---------
        case 'm1'
            m1 = m1_data2struct(m1_UYf, Uf_ivs, Yf_ivs,pctiles);
            clear m1_UYf;
        
        % ----- mLf (Lf values) ---------------------------------------------
        case 'mLf'
            mLf = mLf_data2struct(mLf_data, Cases, nX, pctiles);
            clear mLf_data;
        
        % ----- m3 (DDPC performance) ---------------------------------------
        case 'm3'
            m3 = m3_data2struct(m3_data,cost_types,Cases,pctiles);
            clear m3_data;
        
        % ----- m4 (yf prediction errors) -----------------------------------
        case 'm4'
            m4 = m4_data2struct(m4_data,Cases,nX);
            clear m4_data
    end
end

%% save processed data structures
fndata = 'processed_data.mat';
fprintf('Saving data to %s\n',fndata);
if isfile(fndata)
    save(fndata,OutVars{:},'-append');
else
    save(fndata,OutVars{:});
end

end

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