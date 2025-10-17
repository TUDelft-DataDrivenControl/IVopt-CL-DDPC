function [m0,m1,m2,mLf,m3,m4] = process_dX(data_type,seeds,X_all,opts,plant,OutVars)
arguments
    data_type (1,:) char {mustBeMember(data_type,{'N','Re','p'})}
    seeds double {mustBePositive,mustBeInteger,mustBeMatrix}
    X_all (1,:) double {mustBePositive}
    opts (1,1) struct
    plant ss
    OutVars (1,:) string {mustBeMember(OutVars,{'m0','m1','m2','mLf','m3','m4'})} = {'m0','m1','m2','mLf','m3','m4'};
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
if ~strcmp(data_type,'p')
    Gamma_f = make_ext_obsv(A,C,f);
    Gamma_p = make_ext_obsv(A,C,p);
    Hp = make_blk_tril_toeplitz(A,K,C,eye(ny),p);
    effEpMat = -Gamma_f*(A-K*C)^p*pinv(Gamma_p)*Hp;
end

%% iterate over Re \ N values - initial data processing
subdir1 = pwd; % should be data/sys#/raw/dX/<subdir1>

% get all <subdir2> directories in data/raw/sys#/dX/<subdir1>
subdir2s = dir(subdir1);
isub = [subdir2s(:).isdir]; 
subdir2s = {subdir2s(isub).name};
subdir2s = subdir2s(~ismember(subdir2s,{'.','..','mfiles'}));

% iterate over all Re \ N \ p values
nX = numel(subdir2s);

% find Cases used
seed_mat_files = dir(fullfile(pwd,subdir2s{1},'seed_*.mat'));
Cases = load(seed_mat_files(1).name,'Cases').Cases;
% expect selection of {'iv1','iv2a','iv2b','iv2c','iv3a','iv3c','iv4a','iv4b','iv4c', 'iv5a','iv5b','iv5c','iv6a','iv6c','CLSPC','actLf','TrPred'};
num_Cases = numel(Cases);

%% initializing measures
% ======================== initialize measure 0 (m0) ======================
% -> statistics of Uf & Yf values
% -> structure: m0.<Uf/Yf>.<IVname>.(iX1/2/3...).(mean/median/pctiles)
%     example: m0.Uf.iv1.iX1.mean         (nu, ndiags)
%              m0.Yf.iv2b.iX3.pctiles     (ny, ndiags, num_pctiles)
% -> sliceable cell array containers:
%    m0_Uf_mean     (num_Uf_ivs, nX)  cells of size (nu, ndiags)
%    m0_Yf_mean     (num_Yf_ivs, nX)                (ny, ndiags)
%    m0_Uf_median   (num_Uf_ivs, nX)                (nu, ndiags)
%    m0_Yf_median   (num_Yf_ivs, nX)                (ny, ndiags)
%    m0_Uf_pctiles  (num_Uf_ivs, nX)                (nu, ndiags, num_pctiles)
%    m0_Yf_pctiles  (num_Yf_ivs, nX)                (ny, ndiags, num_pctiles)

% --- IV definitions
Uf_ivs = {'iv1','iv2a','iv2c','iv3c','iv4a','iv4c','iv5a','iv5c','iv6c'}; Uf_ivs = intersect(Uf_ivs,Cases);
Yf_ivs = {'iv2b','iv3a','iv4b','iv5b','iv6a'};                            Yf_ivs = intersect(Yf_ivs,Cases);
num_Uf_ivs = numel(Uf_ivs); % needed for nested for loop inside parfor
num_Yf_ivs = numel(Yf_ivs); % needed for nested for loop inside parfor

% initialize cell arrays
m0_Uf_mean = cell(num_Uf_ivs,nX);    % mean         cell sizes: nu x ndiags
m0_Yf_mean = cell(num_Yf_ivs,nX);    %                          ny x ndiags
m0_Uf_median = m0_Uf_mean;           % median
m0_Yf_median = m0_Yf_mean;
m0_Uf_pctiles = cell(num_Uf_ivs,nX); % percentiles  cell sizes: nu x ndiags x num_pctiles
m0_Yf_pctiles = cell(num_Yf_ivs,nX); %                          ny x ndiags x num_pctiles

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

pctiles = 0:5:100;

% --- Preallocate "sliced" containers ---
m1_Uf_data = zeros(num_Uf_ivs, nX, spX);
m1_Yf_data = zeros(num_Yf_ivs, nX, spX);

% ======================== initialize measure 2 (m2) ======================
% -> identification error (frobenius norm)
% -> structure: m2.<IDerrorType>.<caseName>.data    (nX, spX)
%                                          .mean    (nX, 1)
%                                          .median  (nX, 1)
%                                          .pctiles (nX, num_pctiles)
%      example: m2.Up.iv1.data(iX,ks)
% -> sliceable containers:
%    m2_data(num_IDerrorTypes, num_Cases, nX, spX)

% --- types of identification error
IDerrorTypes = {'Up','Yp','Uf'};
num_IDerrorTypes = numel(IDerrorTypes);

% --- Preallocate sliceable containers ---
m2_data  = zeros(num_IDerrorTypes, num_Cases, nX, spX); % main data
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
m3_data = zeros(numel(cost_types), numel(Cases), nX, spX);   % main data

% ======================== measure 4 (prediction error) ===================
% -> prediction error
[Yf_RelErr_sd,Yf_RelErr_mean] = deal(zeros(nX,ny*f,num_Cases,4));

%% -------------------- loop over Re \ N \ p values -----------------------
% Check and start parallel pool if needed
if ismember('SlurmProfile1',parallel.clusterProfiles) % on cluster?
    if ~isempty(gcp('nocreate'))
        delete(gcp('nocreate'));
    end
    myCluster = parcluster('SlurmProfile1');
    nworkers = 30;
    myCluster.SubmitArguments = SlurmSubmitArgs(['calc_',data_type],15,ntasks=nworkers,cpt=1,GB=3.8);
    parpool(myCluster,nworkers);
else
    if ~isempty(gcp('nocreate'))
        curr_pool = gcp('nocreate');
        if curr_pool.NumWorkers < feature('numcores')
            delete(curr_pool);
            parpool('local',feature('numcores'));
        end
    else
        parpool('local', feature('numcores'));
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
                
                % also calculate matrix describing effect of past noise for m4:
                Gamma_f = make_ext_obsv(A,C,f);
                Gamma_p = make_ext_obsv(A,C,p);
                Hp = make_blk_tril_toeplitz(A,K,C,eye(ny),p);
                effEpMat = -Gamma_f*(A-K*C)^p*pinv(Gamma_p)*Hp;
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
    
    
    % ======================== processing m0, m1, m2, m3 ==================
    fprintf('Processing %s index [%d/%d] (%s = %g) in subdir: %s\n', data_type, iX, nX, data_type, X, subdir2);
    
    % ------------------------ calculations for m0 ------------------------
    if ismember('m0',OutVars)
    [m0_Uf_mean(:,iX),   m0_Yf_mean(:,iX),    m0_Uf_median(:,iX), ...
     m0_Yf_median(:,iX), m0_Uf_pctiles(:,iX), m0_Yf_pctiles(:,iX)] ...
     = process_m0(Uf_ivs,Yf_ivs,iX,nu,ny,f,N,seeds,spX,pctiles);
    end

    % --------------------- calculations for m1,m2,m3 ---------------------
    if any(ismember({'m1','m2','mLf','m3'},OutVars))
    % iterates over noise realizations
    [m1_Uf_data(:, iX, :), m1_Yf_data(:, iX, :),  m2_data(:, :, iX, :), mLf_data{iX},...
     m3_data(:, :, iX, :), Yf_RelErr_sd(iX,:,:,:),Yf_RelErr_mean(iX,:,:,:)] ...
     = process_m123(Uf_ivs,Yf_ivs,Cases,IDerrorTypes,cost_types,iX,nu,ny,p,f,seeds,spX,Hf,effEpMat,OutVars);
    end

    cd(subdir1);
end
if ismember('SlurmProfile1',parallel.clusterProfiles)
    delete(gcp('nocreate')); % close parallel pool
end

%% processing data - m0 (IV statistics)
if ismember('m0',OutVars)
    fprintf("Processing m0 data\n")
    m0 = struct;
    % Uf part
    tic
    for k = 1:numel(Uf_ivs)
        for iX = 1:nX
            iXstr = sprintf('iX%d',iX);
            m0.Uf.(Uf_ivs{k}).(iXstr).mean     = m0_Uf_mean{k,iX};
            m0.Uf.(Uf_ivs{k}).(iXstr).median   = m0_Uf_median{k,iX};
            m0.Uf.(Uf_ivs{k}).(iXstr).pctiles  = m0_Uf_pctiles{k,iX};
        end
    end
    clear m0_Uf_mean m0_Uf_median m0_Uf_pctiles
    toc

    % Yf part
    tic
    for k = 1:numel(Yf_ivs)
        for iX = 1:nX
            iXstr = sprintf('iX%d',iX);
            m0.Yf.(Yf_ivs{k}).(iXstr).mean     = m0_Yf_mean{k,iX};
            m0.Yf.(Yf_ivs{k}).(iXstr).median   = m0_Yf_median{k,iX};
            m0.Yf.(Yf_ivs{k}).(iXstr).pctiles  = m0_Yf_pctiles{k,iX};
        end
    end
    clear m0_Yf_mean m0_Yf_median m0_Yf_pctiles
    toc
else
    m0 = [];
end

%% processing data - m1 (quality of the approximation of the optimal IV)
if ismember('m1',OutVars)
    fprintf("Processing m1 data\n")
    m1 = struct();
    tic
    % Uf part
    for k = 1:numel(Uf_ivs)
        m1.Uf.(Uf_ivs{k}).data     = squeeze( m1_Uf_data(k,:,:) );
        m1.Uf.(Uf_ivs{k}).mean     = mean(    m1.Uf.(Uf_ivs{k}).data, 2);
        m1.Uf.(Uf_ivs{k}).median   = median(  m1.Uf.(Uf_ivs{k}).data, 2);
        m1.Uf.(Uf_ivs{k}).pctiles  = prctile( m1.Uf.(Uf_ivs{k}).data, pctiles, 2);
    end
    clear m1_Uf_data
    toc

    tic
    % Yf part
    for k = 1:numel(Yf_ivs)
        m1.Yf.(Yf_ivs{k}).data     = squeeze( m1_Yf_data(k,:,:) );
        m1.Yf.(Yf_ivs{k}).mean     = mean(    m1.Yf.(Yf_ivs{k}).data, 2);
        m1.Yf.(Yf_ivs{k}).median   = median(  m1.Yf.(Yf_ivs{k}).data, 2);
        m1.Yf.(Yf_ivs{k}).pctiles  = prctile( m1.Yf.(Yf_ivs{k}).data, pctiles, 2);
    end
    clear m1_Yf_data
    toc
else
    m1 = [];
end

%% processing data - m2 (identification error)
if ismember('m2',OutVars)
    fprintf("Processing m2 data\n")
    m2 = struct(); % <- m2_data(kType, kIVn, iX, ks)
    tic
    for kType = 1:numel(IDerrorTypes)
        typeName = IDerrorTypes{kType};
        for kIVn = 1:numel(Cases)
            caseName = Cases{kIVn};
            m2.(typeName).(caseName).data    = squeeze( m2_data(kType, kIVn, :, :) );
            m2.(typeName).(caseName).mean    = mean(    m2.(typeName).(caseName).data, 2);
            m2.(typeName).(caseName).median  = median(  m2.(typeName).(caseName).data, 2);
            m2.(typeName).(caseName).pctiles = prctile( m2.(typeName).(caseName).data, pctiles, 2);
        end
    end
    clear m2_data
    toc
else
    m2 = [];
end

if ismember('mLf',OutVars)
    fprintf("Processing mLf data\n")
    mLf = struct(); % <- mLf_data = nan(nX,num_Cases,spX,ny*f,size(Lf,1),size(Lf,2));
    tic
    for kC = 1:num_Cases
        caseName = Cases{kC};
        for iX = 1:nX
            iXstr = sprintf('iX%d',iX);
            mLf_data_iX_kC = squeeze(mLf_data{iX}(kC, :, :,:));
            if strcmp(caseName,'actLf')
                mLf.(caseName).(iXstr) = squeeze(mLf_data_iX_kC(1,:,:)); % same for all seeds
            else
                mLf.(caseName).(iXstr).data    = mLf_data_iX_kC;
                mLf.(caseName).(iXstr).mean    = squeeze( mean(    mLf_data_iX_kC, 1) );
                mLf.(caseName).(iXstr).median  = squeeze( median(  mLf_data_iX_kC, 1) );
                mLf.(caseName).(iXstr).pctiles = squeeze( prctile( mLf_data_iX_kC, pctiles, 1) );
                mLf.(caseName).(iXstr).std     = squeeze( std( mLf_data_iX_kC, 0 ,1) );
            end
        end
    end
    clear mLf_data
    toc
else
    mLf = [];
end

%% processing data - m3 (DDPC performance)
if ismember('m3',OutVars)
    fprintf("Processing m3 data\n")
    m3 = struct(); % <- m3_data(kType, kIVn, iX, ks)
    tic
    for kType = 1:numel(cost_types)
        costName = cost_types{kType};
        for kIVn = 1:numel(Cases)
            caseName = Cases{kIVn};
            m3.(costName).(caseName).data    = squeeze( m3_data(kType, kIVn, :, :) );
            m3.(costName).(caseName).mean    = mean(    m3.(costName).(caseName).data, 2);
            m3.(costName).(caseName).median  = median(  m3.(costName).(caseName).data, 2);
            m3.(costName).(caseName).pctiles = prctile( m3.(costName).(caseName).data, pctiles, 2);
        end
    end
    clear m3_data
    toc
else
    m3 = [];
end

%% processing data - m4 (yf prediction errors)
% since data can vary greatly in magnitude: use relevative difference
if ismember('m4',OutVars)
    fprintf("Processing m4 data\n")
    m4 = struct;
    tic;
    for kC = 1:numel(Cases)
        caseName = Cases{kC};
        for iX = 1:nX
            iXstr = sprintf('iX%d',iX);
            % (Yf_hat  - Yf)./Yf -> yfhat
            m4.yfhat.(caseName).(iXstr).mean = squeeze(Yf_RelErr_mean(iX,:,kC,1)); % mean over spX seeds
            m4.yfhat.(caseName).(iXstr).std  = squeeze(Yf_RelErr_sd(iX,:,kC,1));   % mean std. dev. per seed over spX seeds

            % (Yf_hatS  - Yf)./Yf -> yfhatS
            m4.yfhatS.(caseName).(iXstr).mean = squeeze(Yf_RelErr_mean(iX,:,kC,2)); % mean over spX seeds
            m4.yfhatS.(caseName).(iXstr).std  = squeeze(Yf_RelErr_sd(iX,:,kC,2));   % mean std. dev. per seed over spX seeds

            % Yf_by_Ep./Yf -> MatEp        ratio of contribution of Ep to Yf
            m4.MatEp.(caseName).(iXstr).mean = squeeze(Yf_RelErr_mean(iX,:,kC,3)); % mean over spX seeds
            m4.MatEp.(caseName).(iXstr).std  = squeeze(Yf_RelErr_sd(iX,:,kC,3));   % mean std. dev. per seed over spX seeds

            % Yf_by_Ef./Yf -> HfEf          = Hf*Ef./Yf
            m4.HfEf.(caseName).(iXstr).mean = squeeze(Yf_RelErr_mean(iX,:,kC,4)); % mean over spX seeds
            m4.HfEf.(caseName).(iXstr).std  = squeeze(Yf_RelErr_sd(iX,:,kC,4));   % mean std. dev. per seed over spX seeds
        end
    end
    toc
else
    m4 = [];
end

%% save processed data
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