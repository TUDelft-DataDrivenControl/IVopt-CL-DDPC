function [m1_UYf_iX, mLf_data_iX, m4_data_iX] ...
    = process_m_all(Uf_ivs,Yf_ivs,Cases,iX,nu,ny,p,f,seeds,spX,Hf,OutVars)
%% ======================== initialize data containers ======================

num_OutVars = numel(OutVars);
% initializing max counters needed for nested for loop inside parfor
num_Uf_ivs       = numel(Uf_ivs); 
num_Yf_ivs       = numel(Yf_ivs);
num_Cases        = numel(Cases);
num_IVs          = num_Uf_ivs + num_Yf_ivs;

m1_UYf_iX     = zeros(num_IVs, spX);
cols_Lf = (ny+nu)*p+ny*f;
mLf_data_iX   = zeros(num_Cases,spX,ny*f,cols_Lf);

% measure 4: future output relative errors: actual, predictions & contributions (by Ep & Ef)
% Combine the 4 error components into one 4-D array: (ny*f) x num_Cases x spX x 3
m4_data_iX = zeros(ny*f, num_Cases, spX, 3); % (:,:,:,i): i=1 -> std. dev, i=2 -> mean, i=3 -> rms

% determine if running on cluster
if ismember('SlurmProfile1',parallel.clusterProfiles) % running on cluster?
    onCluster = true;
else
    onCluster = false;
end

% waitbar setup
if ~onCluster
    % Initialize progress tracking
    w = waitbar(0,sprintf('Progress: (%d/%d) -> %.1f%%',0,spX,0));
    D = parallel.pool.DataQueue;
    afterEach(D,@parforWaitbar);
    parforWaitbar(w, spX); % Initialize waitbar function
else
    D = []; % needed to prevent error in parfor
end

%% ---------------- Parallel iteration over seeds -------------------------
m_tictoc = tic;
txt_iters = sprintf('\tIterating over seeds to calculate measures.');
fprintf('%s',txt_iters);
parfor ks = 1:spX
    seed = seeds(ks,iX);
    fndata = sprintf('seed_%d.mat',seed);
    [Z,u0,y0,e0,e1,u_cl,y_cl,Lf] = load_seedmat(fndata,OutVars);

    % iterating over output variables to calculate
    for kOutVar = 1:num_OutVars
        OutVar = OutVars{kOutVar};
    switch OutVar
    % =================================================================
    % m1) how well IV approximates optimal one
        case 'm1'
            m1_UYf_iX(:, ks) = process_m1_seed(Uf_ivs, Yf_ivs, Z, nu, ny, f);

    % =================================================================
    % mLf) average identified Lf
        case 'mLf'
            mLf_data_iX(:,ks,:,:) = process_mLf_seed(Cases,ny,f,cols_Lf,Lf,ks);

    % =================================================================
    % m4) prediction error
        case 'm4'
            m4_data_iX(:,:,ks,:) = process_m4_seed(e0,e1,u0,y0,u_cl,y_cl,Lf,Hf,Cases,p,f);
    end
    end

    % =================================================================
    if ~onCluster; send(D, []); end % update waitbar
end % of parfor

%% processing of seeds
% reformat data for m4
if ismember('m4',OutVars)
    % mean over seeds (3rd dim) -> squeeze to ny*f x num_Cases x 4 x 2
    m4_data_iX = squeeze(mean(m4_data_iX, 3));
else
    m4_data_iX = zeros(ny*f, num_Cases, 3);
end

m_time = toc(m_tictoc);
if ~onCluster; close(w); end
fprintf('\tFinished in %.2f seconds\n', m_time);

end

%% Helper functions
function [Z,u0,y0,e0,e1,u_cl,y_cl,Lf] = load_seedmat(fnpath,OutVars)
    [Z,u0,y0,e0,e1,u_cl,y_cl,Lf] = deal([]); % will be filled as needed
    % Only load necessary variables to improve memory usage and performance
    loadVars = {};
    for kOutVar = 1:numel(OutVars)
        switch OutVars{kOutVar}
            case 'm1'
                loadVars = [loadVars, 'Z'];
            case 'm4'
                loadVars = [loadVars, 'u0', 'y0', 'e0', 'e1', 'u_cl', 'y_cl', 'Lf'];
            case 'mLf'
                loadVars = [loadVars, 'Lf'];
        end
    end
    loadVars = unique(loadVars); % avoid duplicates (see 'Lf' for m4 and mLf)
    s = load(fnpath,loadVars{:});
    for kVar = 1:numel(loadVars)
        varName = loadVars{kVar};
        switch varName
            case 'Z'
                Z = s.Z;
            case 'u0'
                u0 = s.u0;
            case 'y0'
                y0 = s.y0;
            case 'e0'
                e0 = s.e0;
            case 'e1'
                e1 = s.e1;
            case 'u_cl'
                u_cl = s.u_cl;
            case 'y_cl'
                y_cl = s.y_cl;
            case 'Lf'
                Lf = s.Lf;
        end
    end
end

function parforWaitbar(waitbarHandle,iterations)
    persistent count h Nit
    
    if nargin == 2
        % Initialize
        
        count = 0;
        h = waitbarHandle;
        Nit = iterations;
    else
        % Update the waitbar
        
        % Check whether the handle is a reference to a deleted object
        if isvalid(h)
            count = count + 1;
            waitbar(count / Nit,h,sprintf('Progress: (%d/%d) -> %.1f%%',count,Nit,count/Nit*100));
        end
    end
end