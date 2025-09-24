%% iterate over N values - initial data processing
% get all <subdir2> directories
subdir2s = dir(pwd);
isub = [subdir2s(:).isdir]; 
subdir2s = {subdir2s(isub).name};
subdir2s = subdir2s(~ismember(subdir2s,{'.','..','mfiles'}));

% iterate over all N values
nN = numel(subdir2s);
nX = nN; % for scripts performing calculations (agnostic to varying Re or N)

%% initializing measures
% ======================== initialize measure 0 (m0) ======================
% Uf & Yf values

% --- IV definitions
Uf_ivs = {'iv1','iv2a','iv2c','iv3c','iv4a','iv4c','iv5a','iv5c','iv6c'};
Yf_ivs = {'iv2b','iv3a','iv4b','iv5b','iv6a'};
num_Uf_ivs = numel(Uf_ivs); % needed for nested for loop inside parfor
num_Yf_ivs = numel(Yf_ivs); % needed for nested for loop inside parfor

% initialize cell arrays
m0_Uf_mean = cell(num_Uf_ivs,nX);                % mean         cell sizes: nu x ndiags
m0_Yf_mean = cell(num_Yf_ivs,nX);                %                          ny x ndiags
m0_Uf_median = m0_Uf_mean;                       % median
m0_Yf_median = m0_Yf_mean;
m0_Uf_pctiles = cell(num_Uf_ivs,nX);             % percentiles  cell sizes: nu x ndiags x num_pctiles
m0_Yf_pctiles = cell(num_Yf_ivs,nX);             %                          ny x ndiags x num_pctiles

% ======================== initialize measure 1 (m1) ======================
% -> m1: how well IV approximates optimal IV

pctiles = 0:5:100;
num_pctiles  = numel(pctiles);

% --- Preallocate "sliced" containers ---
m1_Uf_data     = zeros(num_Uf_ivs, nX, spX);
m1_Yf_data     = zeros(num_Yf_ivs, nX, spX);

% other
iyf = nu*f + (1:ny*f);

% ======================== initialize measure 2 (m2) ======================
% -> identification error (frobenius norm)

% --- IV/case names
Cases = {'iv1','iv2a','iv2b','iv2c','iv3a','iv3c','iv4a','iv4b','iv4c', ...
         'iv5a','iv5b','iv5c','iv6a','iv6c','CLSPC','actLf'};
num_Cases = numel(Cases);


IDerrorTypes = {'Up','Yp','Uf'};  % types of identification error

% --- Preallocate sliceable containers ---
m2_data    = zeros(numel(IDerrorTypes), numel(Cases), nX, spX); % main data

num_IDerrorTypes = numel(IDerrorTypes);

% ======================== initialize measure 3 (m3) ======================
% -> DDPC performance

% --- cost types and cases
cost_types = {'cost_u','cost_y','cost_tot'};
num_cost_types = numel(cost_types);

% --- Preallocate sliceable containers ---
m3_data    = zeros(numel(cost_types), numel(Cases), nX, spX);   % main data

%% --------------------------- loop over N values -------------------------
for iN = 1:nN
    N = N_all(iN);

    iX = iN; % for scripts performing calculations

    % choose subdir2 corresponding with N value
    subdir2 = choose_subdir_by_number(subdir2s, N);
    cd(subdir2);                                % navigate into <subdir2>
    
    % load file with settings <subdir2>
    settingsFile = find_settingsFile();
    load(settingsFile);
    
    % ======================== calculations for m0 ========================
    num_diags = f+N-1;
               
    % ---------------------------- for Uf ---------------------------------
    m0_Uf_tictoc = tic;
    parfor kIVu = 1:num_Uf_ivs
        iv_name = Uf_ivs{kIVu};

        % initialize m0_Uf_...
        m0_Uf_mean2 = zeros(nu,num_diags);                % mean
        m0_Uf_median2 = m0_Uf_mean2;                      % median
        m0_Uf_pctiles2 = zeros(nu,num_diags,num_pctiles); % percentiles

        % get Uf values in all seeds
        m0_Uf_data2 = zeros(nu*f,N,spX); % initialize
        for ks = 1:spN
            seed = seeds(iN,ks);
            fndata = sprintf('seed_%d.mat',seed);
            Z = load(fndata,'Z').Z;
            
            Uf_iv = Z.([iv_name,'_']);
            m0_Uf_data2(:,:,ks) = Uf_iv;
        end
        
        % get all values belonging to anti-diagonals over seeds
        ij_adiags = get_subind_diags(f,N,nr=nu,anti=true); % i & j indices for 2D case
        for kd = 1:num_diags
            rows = ij_adiags{kd}(:,1);  % row indices for 2D case
            cols = ij_adiags{kd}(:,2);  % col indices for 2D case

            % get linear index for 3D case
            rows2 = repmat(rows,spX,1); % extend (spX seeds)
            cols2 = repmat(cols,spX,1);
            i3s = kron( (1:spX).', ones(numel(cols),1) ); % 3rd dim indices
            idxlin = sub2ind(size(m0_Uf_data2),rows2,cols2,i3s); % convert to linear indices
            Uf_sel = m0_Uf_data2(idxlin);
            Uf_sel = reshape(Uf_sel,nu,[]);

        % calculate mean, median, percentiles
            m0_Uf_mean2(:,kd)      = mean(Uf_sel,2);
            m0_Uf_median2(:,kd)    = median(Uf_sel,2);
            m0_Uf_pctiles2(:,kd,:) = prctile(Uf_sel,pctiles,2);
        end
        
        % assign cells in data arrays
        m0_Uf_mean{kIVu,iX}    = m0_Uf_mean2;
        m0_Uf_median{kIVu,iX}  = m0_Uf_median2;
        m0_Uf_pctiles{kIVu,iX} = m0_Uf_pctiles2;
    end
    m0_Uf_time = toc(m0_Uf_tictoc);
    fprintf("m0 for Uf finished in %.2f seconds\n",m0_Uf_time);

    % ---------------------------- for Yf ---------------------------------
    m0_Yf_tictoc = tic;
    parfor kIVy = 1:num_Yf_ivs
        iv_name = Yf_ivs{kIVy};
        
        % initialize m0_Yf_...
        m0_Yf_mean2 = zeros(ny,num_diags);                % mean
        m0_Yf_median2 = m0_Yf_mean2;                      % median
        m0_Yf_pctiles2 = zeros(ny,num_diags,num_pctiles); % percentiles

        % get Yf values in all seeds
        m0_Yf_data2 = zeros(ny*f,N,spX); % initialize
        for ks = 1:spN
            seed = seeds(iN,ks);
            fndata = sprintf('seed_%d.mat',seed);
            Z = load(fndata,'Z').Z;
            
            switch iv_name
                case 'iv3a' % IV_Theta: only possible because ny = nlcf (see get_Z.m)
                    Yf_iv = Z.iv3a_(1:ny*f,:);

                case 'iv6a' % Rf_yr0 (future references)
                    Yf_iv = Z.iv6a_;

                otherwise
                    Yf_iv = Z.([iv_name,'_'])(iyf,:);
            end
            m0_Yf_data2(:,:,ks) = Yf_iv;
        end
        
        % get all values belonging to anti-diagonals over seeds
        ij_adiags = get_subind_diags(f,N,nr=ny,anti=true); % i & j indices for 2D case
        for kd = 1:num_diags
            rows = ij_adiags{kd}(:,1);  % row indices for 2D case
            cols = ij_adiags{kd}(:,2);  % col indices for 2D case

            % get linear index for 3D case
            rows2 = repmat(rows,spX,1); % extend (spX seeds)
            cols2 = repmat(cols,spX,1);
            i3s = kron( (1:spX).', ones(numel(cols),1) ); % 3rd dim indices
            idxlin = sub2ind(size(m0_Yf_data2),rows2,cols2,i3s); % convert to linear indices
            Yf_sel = m0_Yf_data2(idxlin);
            Yf_sel = reshape(Yf_sel,ny,[]);

        % calculate mean, median, percentiles
            m0_Yf_mean2(:,kd)      = mean(Yf_sel,2);
            m0_Yf_median2(:,kd)    = median(Yf_sel,2);
            m0_Yf_pctiles2(:,kd,:) = prctile(Yf_sel,pctiles,2);
        end

        % assign cells in data arrays
        m0_Yf_mean{kIVy,iX}    = m0_Yf_mean2;
        m0_Yf_median{kIVy,iX}  = m0_Yf_median2;
        m0_Yf_pctiles{kIVy,iX} = m0_Yf_pctiles2;
    end
    m0_Yf_time = toc(m0_Yf_tictoc);
    fprintf("m0 for Yf finished in %.2f seconds\n",m0_Yf_time);

    % ===================== calculations for m1,m2,m3 =====================
    % iterate over noise realizations
    tic
    parfor ks = 1:spN
        seed = seeds(iN,ks);
        fndata = sprintf('seed_%d.mat',seed);
        [Cases,Cz,FroIDerror,Lf,Tcl,Z,cost_tot,cost_u,cost_u1,cost_u2,cost_y,...
          e0,e1,opts,u0,u_cl,u_iv,xcl0,y0,y_cl,y_iv] = load_seedmat(fndata);
        
        % =================================================================
        % m1) how well IV approximates optimal one
        %   -> calculates m1_Uf_data(kIVu, iX, ks) and m1_Yf_data(kIVy, iX, ks)

        % ---- Uf IVs ----
        for kIVu = 1:num_Uf_ivs
            iv_name = Uf_ivs{kIVu};
            
            switch iv_name
                case 'iv2a'
                    % leave zero since this is the optimal IV for Uf
                otherwise
                    Uf_iv = Z.([iv_name,'_']);
                    m1_Uf_data(kIVu, iX, ks) = norm(Uf_iv - Z.iv2a_, 'fro');
            end
        end

        % ---- Yf IVs ----
        for kIVy = 1:num_Yf_ivs
            iv_name = Yf_ivs{kIVy};

            switch iv_name
                case 'iv2b'
                    % leave zero since this contains the optimal IV for Yf
                    Yf_iv = 0; % simply to suppress warning
                    calc_norm = false;

                case 'iv3a' % IV_Theta: only possible because ny = nlcf (see get_Z.m)
                    Yf_iv = Z.iv3a_(1:ny*f,:);
                    calc_norm = true;

                case 'iv6a' % Rf_yr0 (future references)
                    Yf_iv = Z.iv6a_;
                    calc_norm = true;

                otherwise
                    Yf_iv = Z.([iv_name,'_'])(iyf,:);
                    calc_norm = true;
            end

            if calc_norm
                m1_Yf_data(kIVy, iX, ks) = norm(Yf_iv - Z.iv2b_(iyf,:), 'fro');
            end
        end

        % =================================================================
        % m2) identification error (frobenius norm)
        
        for kType = 1:num_IDerrorTypes         %  1   2   3
            IDerrorType = IDerrorTypes{kType}; % Up, Yp, Uf

            for kIVn = 1:num_Cases
                IVn = Cases{kIVn}; % iv1, iv2a, CLSPC, etc.
                m2_data(kType, kIVn, iX, ks) = FroIDerror.(IVn).(IDerrorType);
            end
        end

        % =================================================================
        % m3) DDPC performance
        
        for kType = 1:num_cost_types
            cost_type = cost_types{kType};

            % pick the right cost struct
            switch cost_type
                case 'cost_u'
                    cost = cost_u;
                case 'cost_y'
                    cost = cost_y;
                otherwise %'cost_tot'
                    cost = cost_tot;
            end

            for kIVn = 1:num_Cases
                IVn = Cases{kIVn};
                m3_data(kType, kIVn, iX, ks) = cost.(IVn);
            end
        end
        % =================================================================
        
    end % of parfor
    toc

    cd(subdir1);
end

%% processing data - m0 (IV statistics)
fprintf("Processing m0 data\n")

m0 = struct;
% Uf part
for k = 1:numel(Uf_ivs)
    for iX = 1:nX
        iXstr = sprintf('iX%d',iX);
        m0.Uf.(Uf_ivs{k}).(iXstr).mean     = m0_Uf_mean{k,iX};
        m0.Uf.(Uf_ivs{k}).(iXstr).median   = m0_Uf_median{k,iX};
        m0.Uf.(Uf_ivs{k}).(iXstr).pctiles  = m0_Uf_pctiles{k,iX};
    end
end
% clear m0_Uf_mean m0_Uf_median m0_Uf_pctiles

% Yf part
for k = 1:numel(Yf_ivs)
    for iX = 1:nX
        iXstr = sprintf('iX%d',iX);
        m0.Yf.(Yf_ivs{k}).(iXstr).mean     = m0_Yf_mean{k,iX};
        m0.Yf.(Yf_ivs{k}).(iXstr).median   = m0_Yf_median{k,iX};
        m0.Yf.(Yf_ivs{k}).(iXstr).pctiles  = m0_Yf_pctiles{k,iX};
    end
end
% clear m0_Yf_mean m0_Yf_median m0_Yf_pctiles

%% processing data - m1 (quality of the approximation of the optimal IV)
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

%% processing data - m2 (identification error)
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

%% processing data - m3 (DDPC performance)
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

%% save processed data
fndata = 'processed_data.mat';
fprintf('Saving data to %s\n',fndata);
save(fndata,'m0','m1','m2','m3');

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

function [Cases,Cz,FroIDerror,Lf,Tcl,Z,cost_tot,cost_u,cost_u1,cost_u2,cost_y,...
          e0,e1,opts,u0,u_cl,u_iv,xcl0,y0,y_cl,y_iv] = load_seedmat(fnpath)
    load(fnpath);
end