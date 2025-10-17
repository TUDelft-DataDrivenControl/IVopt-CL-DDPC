function [m1_Uf_data_iX, m1_Yf_data_iX, m2_data_iX,mLf_data_iX, m3_data_iX,...
    Yf_RelErr_sd_iX,Yf_RelErr_mean_iX] ...
    = process_m123(Uf_ivs,Yf_ivs,Cases,IDerrorTypes,cost_types,iX,nu,ny,p,f,seeds,spX,Hf,effEpMat,OutVars)
num_OutVars = numel(OutVars);
% initializing max counters needed for nested for loop inside parfor
num_Uf_ivs       = numel(Uf_ivs); 
num_Yf_ivs       = numel(Yf_ivs); 
num_Cases        = numel(Cases);
num_IDerrorTypes = numel(IDerrorTypes);
num_cost_types   = numel(cost_types);

m1_Uf_data_iX = zeros(num_Uf_ivs,spX);
m1_Yf_data_iX = zeros(num_Yf_ivs,spX);
m2_data_iX    = zeros(num_IDerrorTypes,num_Cases,spX);
m3_data_iX    = zeros(num_cost_types,num_Cases,spX);
cols_Lf = (ny+nu)*p+ny*f;
mLf_data_iX   = zeros(num_Cases,spX,ny*f,cols_Lf);

% measure 4: future output errors: actual, predictions & error contributions (by Ep & Ef)
% Combine the 4 error components into one 4-D array: (ny*f) x num_Cases x spX x 4
Yf_RelErr_sd_iX   = zeros(ny*f, num_Cases, spX, 4);
Yf_RelErr_mean_iX = zeros(ny*f, num_Cases, spX, 4);

if ismember('SlurmProfile1',parallel.clusterProfiles) % running on cluster?
    onCluster = true;
else
    onCluster = false;
end

if ~onCluster
    % Initialize progress tracking
    w = waitbar(0,sprintf('Progress: (%d/%d) -> %.1f%%',0,spX,0));
    D = parallel.pool.DataQueue;
    afterEach(D,@parforWaitbar);
    parforWaitbar(w, spX); % Initialize waitbar function
else
    D = []; % needed to prevent error in parfor
end

m123_tictoc = tic;
txt_iters = sprintf('\tm1, m2, m3: Iterating over seeds (see waitbar).');
fprintf('%s',txt_iters);
parfor ks = 1:spX
    seed = seeds(ks,iX);
    fndata = sprintf('seed_%d.mat',seed);
    [FroIDerror,Z,cost_tot,cost_u,cost_y,u0,y0,e0,e1,u_cl,y_cl,Lf] = load_seedmat(fndata,OutVars);

    % iterating over output variables to calculate
    for kOutVar = 1:num_OutVars
        OutVar = OutVars{kOutVar};
    switch OutVar
    
    % =================================================================
    % m1) how well IV approximates optimal one
        case 'm1'
            [m1_Uf_data_iX(:,ks), m1_Yf_data_iX(:,ks)] = process_m1_seed(Uf_ivs, Yf_ivs, Z, nu, ny, f);

    % =================================================================
    % m2) identification error (frobenius norm)
        case 'm2'
            m2_data_iX(:,:,ks) = process_m2_seed(FroIDerror, Cases, IDerrorTypes);

    % mLf) average identified Lf
        case 'mLf'
            mLf_data_iX(:,ks,:,:) = process_mLf_seed(Cases,ny,f,cols_Lf,Lf,ks);

    % =================================================================
    % m3) DDPC performance
        case 'm3'
            m3_data_iX(:,:,ks) = process_m3_seed(cost_tot,cost_u,cost_y,Cases,cost_types);

    % =================================================================
    % m4) prediction error
        case 'm4'
            [Yf_RelErr_sd_iX(:,:,ks,:), Yf_RelErr_mean_iX(:,:,ks,:)] = process_m4_seed(e0,e1,u0,y0,u_cl,y_cl,Lf,effEpMat,Hf,Cases,p,f);

    end
    end

    % =================================================================
    if ~onCluster; send(D, []); end % update waitbar
end % of parfor

% reformat data for m4
if ismember('m4',OutVars)
    % mean over seeds (3rd dim). Result: ny*f x num_Cases x 1 x 4 -> squeeze to ny*f x num_Cases x 4
    Yf_RelErr_mean_iX = squeeze(mean(Yf_RelErr_mean_iX, 3));
    Yf_RelErr_sd_iX   = squeeze(mean(Yf_RelErr_sd_iX,   3));
end

m123_time = toc(m123_tictoc);
if ~onCluster; close(w); end
fprintf([repmat('\b',1,numel(txt_iters)),'\tm1, m2, m3: Iterating over seeds\t\tFinished in %.2f seconds\n'], m123_time);

end

%% Helper functions
function [FroIDerror,Z,cost_tot,cost_u,cost_y,u0,y0,e0,e1,u_cl,y_cl,Lf] = load_seedmat(fnpath,OutVars)
    [FroIDerror,Z,cost_tot,cost_u,cost_y,u0,y0,e0,e1,u_cl,y_cl,Lf] = deal([]); % will be filled as needed
    % Only load necessary variables to improve memory usage and performance
    loadVars = {};
    for kOutVar = 1:numel(OutVars)
        switch OutVars{kOutVar}
            case 'm1'
                loadVars = [loadVars, 'Z'];
            case 'm2'
                loadVars = [loadVars, 'FroIDerror'];
            case 'm3'
                loadVars = [loadVars, 'cost_tot', 'cost_u', 'cost_y'];
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
            case 'FroIDerror'
                FroIDerror = s.FroIDerror;
            case 'Z'
                Z = s.Z;
            case 'cost_tot'
                cost_tot = s.cost_tot;
            case 'cost_u'
                cost_u = s.cost_u;
            case 'cost_y'
                cost_y = s.cost_y;
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
    persistent count h N
    
    if nargin == 2
        % Initialize
        
        count = 0;
        h = waitbarHandle;
        N = iterations;
    else
        % Update the waitbar
        
        % Check whether the handle is a reference to a deleted object
        if isvalid(h)
            count = count + 1;
            waitbar(count / N,h,sprintf('Progress: (%d/%d) -> %.1f%%',count,N,count/N*100));
        end
    end
end