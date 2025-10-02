%% iterate over Re values - initial data processing
% get all <subdir2> directories in data/raw/dRe/<subdir1>
subdir2s = dir(pwd);
isub = [subdir2s(:).isdir]; 
subdir2s = {subdir2s(isub).name};
subdir2s = subdir2s(~ismember(subdir2s,{'.','..','mfiles'}));

% iterate over all Re values
nX = numel(subdir2s);

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
Uf_ivs = {'iv1','iv2a','iv2c','iv3c','iv4a','iv4c','iv5a','iv5c','iv6c'};
Yf_ivs = {'iv2b','iv3a','iv4b','iv5b','iv6a'};
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
num_pctiles  = numel(pctiles);

% --- Preallocate "sliced" containers ---
m1_Uf_data = zeros(num_Uf_ivs, nX, spX);
m1_Yf_data = zeros(num_Yf_ivs, nX, spX);

% other
iyf = nu*f + (1:ny*f);

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

% --- IV/case names
Cases = {'iv1','iv2a','iv2b','iv2c','iv3a','iv3c','iv4a','iv4b','iv4c', ...
         'iv5a','iv5b','iv5c','iv6a','iv6c','CLSPC','actLf'};
num_Cases = numel(Cases);

% --- Preallocate sliceable containers ---
m2_data = zeros(num_IDerrorTypes, num_Cases, nX, spX); % main data

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
nCostTypes = numel(cost_types);

% --- Preallocate sliceable containers ---
m3_data    = zeros(nCostTypes, numel(Cases), nX, spX);

%% --------------------------- loop over Re values ------------------------
for iX = 1:nX
    Re = Re_all(iX);

    % choose subdir2 corresponding with Re value
    subdir2 = choose_subdir_by_iX(subdir2s, iX);
    cd(subdir2);                                % navigate into <subdir2>
    
    % load file with settings <subdir2>
    % settingsFile = find_settingsFile();
    % load(settingsFile);
    
    % ======================== processing m0, m1, m2, m3 ==================
    fprintf('Processing Re index [%d/%d] (Re = %g) in subdir: %s\n', iX, nX, Re, subdir2);
    
    % ------------------------ calculations for m0 ------------------------
    [m0_Uf_mean(:,iX), m0_Yf_mean(:,iX), m0_Uf_median(:,iX), m0_Yf_median(:,iX), m0_Uf_pctiles(:,iX), m0_Yf_pctiles(:,iX)] = process_m0(Uf_ivs,Yf_ivs,iX,nu,ny,f,N,seeds,spX,pctiles);

    % --------------------- calculations for m1,m2,m3 ---------------------
    % iterates over noise realizations
    [m1_Uf_data(:, iX, :), m1_Yf_data(:, iX, :), m2_data(:, :, iX, :), m3_data(:, :, iX, :)] = process_m123(Uf_ivs,Yf_ivs,Cases,IDerrorTypes,cost_types,iX,nu,ny,f,seeds,spX);

    cd(subdir1);
end

fprintf("Processing m0 data\n")
m0 = struct;
for k = 1:numel(Uf_ivs)
    for iX = 1:nX
        iXstr = sprintf('iX%d',iX);
        m0.Uf.(Uf_ivs{k}).(iXstr).mean     = m0_Uf_mean{k,iX};
        m0.Uf.(Uf_ivs{k}).(iXstr).median   = m0_Uf_median{k,iX};
        m0.Uf.(Uf_ivs{k}).(iXstr).pctiles  = m0_Uf_pctiles{k,iX};
    end
end
for k = 1:numel(Yf_ivs)
    for iX = 1:nX
        iXstr = sprintf('iX%d',iX);
        m0.Yf.(Yf_ivs{k}).(iXstr).mean     = m0_Yf_mean{k,iX};
        m0.Yf.(Yf_ivs{k}).(iXstr).median   = m0_Yf_median{k,iX};
        m0.Yf.(Yf_ivs{k}).(iXstr).pctiles  = m0_Yf_pctiles{k,iX};
    end
end

fprintf("Processing m1 data\n")
m1 = struct();
for k = 1:numel(Uf_ivs)
    m1.Uf.(Uf_ivs{k}).data     = squeeze( m1_Uf_data(k,:,:) );
    m1.Uf.(Uf_ivs{k}).mean     = mean(    m1.Uf.(Uf_ivs{k}).data, 2);
    m1.Uf.(Uf_ivs{k}).median   = median(  m1.Uf.(Uf_ivs{k}).data, 2);
    m1.Uf.(Uf_ivs{k}).pctiles  = prctile( m1.Uf.(Uf_ivs{k}).data, pctiles, 2);
end
for k = 1:numel(Yf_ivs)
    m1.Yf.(Yf_ivs{k}).data     = squeeze( m1_Yf_data(k,:,:) );
    m1.Yf.(Yf_ivs{k}).mean     = mean(    m1.Yf.(Yf_ivs{k}).data, 2);
    m1.Yf.(Yf_ivs{k}).median   = median(  m1.Yf.(Yf_ivs{k}).data, 2);
    m1.Yf.(Yf_ivs{k}).pctiles  = prctile( m1.Yf.(Yf_ivs{k}).data, pctiles, 2);
end

fprintf("Processing m2 data\n")
m2 = struct();
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

fprintf("Processing m3 data\n")
m3 = struct();
for kType = 1:nCostTypes
    costName = cost_types{kType};
    for kIVn = 1:numel(Cases)
        caseName = Cases{kIVn};
        m3.(costName).(caseName).data    = squeeze( m3_data(kType, kIVn, :, :) );
        m3.(costName).(caseName).mean    = mean(    m3.(costName).(caseName).data, 2);
        m3.(costName).(caseName).median  = median(  m3.(costName).(caseName).data, 2);
        m3.(costName).(caseName).pctiles = prctile( m3.(costName).(caseName).data, pctiles, 2);
    end
end

fndata = 'processed_data.mat';
fprintf('Saving data to %s\n',fndata);
save(fndata,'m0','m1','m2','m3');


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
