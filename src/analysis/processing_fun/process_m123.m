function [m1_Uf_data_iX, m1_Yf_data_iX, m2_data_iX,mLf_data_iX, m3_data_iX,...
    Yf_RelErr_sd_iX,Yf_RelErr_mean_iX] ...
    = process_m123(Uf_ivs,Yf_ivs,Cases,IDerrorTypes,cost_types,iX,nu,ny,p,f,seeds,spX,Hf,effEpMat)

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
[Yf_RelErr_sd_iX1,Yf_RelErr_mean_iX1] = deal(zeros(ny*f,num_Cases,spX)); % (Yf_hat  - Yf)./Yf
[Yf_RelErr_sd_iX2,Yf_RelErr_mean_iX2] = deal(zeros(ny*f,num_Cases,spX)); % (Yf_hatS - Yf)./Yf
[Yf_RelErr_sd_iX3,Yf_RelErr_mean_iX3] = deal(zeros(ny*f,num_Cases,spX)); % Yf_by_Ep./Yf
[Yf_RelErr_sd_iX4,Yf_RelErr_mean_iX4] = deal(zeros(ny*f,num_Cases,spX)); % Yf_by_Ef./Yf

iyf = nu*f + (1:ny*f); % row indices of Yf in Z

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
    [FroIDerror,Z,cost_tot,cost_u,cost_y,u0,y0,e0,e1,u_cl,y_cl,Lf] = load_seedmat(fndata);
    
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
                m1_Uf_data_iX(kIVu, ks) = norm(Uf_iv - Z.iv2a_, 'fro');
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
            m1_Yf_data_iX(kIVy, ks) = norm(Yf_iv - Z.iv2b_(iyf,:), 'fro');
        end
    end

    % =================================================================
    % m2) identification error (frobenius norm)
    for kType = 1:num_IDerrorTypes         %  1   2   3
        IDerrorType = IDerrorTypes{kType}; % Up, Yp, Uf

        for kIVn = 1:num_Cases
            IVn = Cases{kIVn}; % iv1, iv2a, CLSPC, etc.
            m2_data_iX(kType, kIVn, ks) = FroIDerror.(IVn).(IDerrorType);
        end
    end

    % mLf) average identified Lf
    mLf_data_iX_ks = zeros(num_Cases,ny*f,cols_Lf);
    for kC = 1:num_Cases
        CaseName = Cases{kC};
        if ~(strcmp(CaseName,'actLf') && ks > 1)
            mLf_data_iX_ks(kC,:,:) = Lf.(CaseName);
        end
    end
    mLf_data_iX(:,ks,:,:) = mLf_data_iX_ks;


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
            m3_data_iX(kType, kIVn, ks) = cost.(IVn);
        end
    end

    % =================================================================
    % m4) prediction error
    % Yf_RelErr_sd : (nX, ny*f, num_Cases, spX, 4)

    e_all = [e0(:,end-p+1:end) e1];
    [~, Ep, Ef] = make_Hankel(e_all, p, f);
    Yf_by_Ep = effEpMat*Ep; % past noise contribution to Yf
    Yf_by_Ef = Hf*Ef;       % future noise contribution to Yf
    for kC = 1:num_Cases
        caseName = Cases{kC};

        % create Hankel matrices
        u_all = [u0(:,end-p+1:end) u_cl.(caseName)];
        [~, Up, Uf] = make_Hankel(u_all, p, f);
        y_all = [y0(:,end-p+1:end) y_cl.(caseName)];
        [~, Yp, Yf] = make_Hankel(y_all, p, f);

        % calculate output predictions
        Yf_hat = Lf.(caseName)*[Up;Yp;Uf];     % prediction w/ Lf estimate
        if ~strcmp(caseName,'actLf')
            Yf_hatS = Lf.('actLf')*[Up;Yp;Uf]; % prediction w/ actual Lf
        else
            Yf_hatS = Yf_hat;
        end
        [Yf_RelErr_sd_iX1(:,kC,ks),Yf_RelErr_mean_iX1(:,kC,ks)] = std( (Yf_hat  - Yf)./Yf, 0, 2); % std. dev & mean over Ncl-f+1 data points
        [Yf_RelErr_sd_iX2(:,kC,ks),Yf_RelErr_mean_iX2(:,kC,ks)] = std( (Yf_hatS - Yf)./Yf, 0, 2);
        [Yf_RelErr_sd_iX3(:,kC,ks),Yf_RelErr_mean_iX3(:,kC,ks)] = std( Yf_by_Ep./Yf, 0, 2);
        [Yf_RelErr_sd_iX4(:,kC,ks),Yf_RelErr_mean_iX4(:,kC,ks)] = std( Yf_by_Ef./Yf, 0, 2);
        % NB: for last dim, indxs 2, 3 & 4 should sum to zero for the mean
    end

    % =================================================================
    if ~onCluster; send(D, []); end % update waitbar
end % of parfor

% reformat data for m4
Yf_RelErr_mean_iX = squeeze(mean(cat(4,Yf_RelErr_mean_iX1,Yf_RelErr_mean_iX2,Yf_RelErr_mean_iX3,Yf_RelErr_mean_iX4),3)); % ny*f, num_Cases,1,4
clear Yf_RelErr_mean_iX1 Yf_RelErr_mean_iX2 Yf_RelErr_mean_iX3 Yf_RelErr_mean_iX4
Yf_RelErr_sd_iX  = squeeze(mean(cat(4,Yf_RelErr_sd_iX1, Yf_RelErr_sd_iX2, Yf_RelErr_sd_iX3, Yf_RelErr_sd_iX4) ,3));
clear Yf_RelErr_std_iX1 Yf_RelErr_std_iX2 Yf_RelErr_std_iX3 Yf_RelErr_std_iX4

m123_time = toc(m123_tictoc);
if ~onCluster; close(w); end
fprintf([repmat('\b',1,numel(txt_iters)),'\tm1, m2, m3: Iterating over seeds\t\tFinished in %.2f seconds\n'], m123_time);

end

%% Helper functions
function [FroIDerror,Z,cost_tot,cost_u,cost_y,u0,y0,e0,e1,u_cl,y_cl,Lf] = load_seedmat(fnpath)
    % Only load necessary variables to improve memory usage and performance
    s = load(fnpath,'FroIDerror','Z','cost_tot','cost_u','cost_y','u0','y0','e0','e1','u_cl','y_cl','Lf');
    FroIDerror = s.FroIDerror;
    Z = s.Z;
    cost_tot = s.cost_tot;
    cost_u = s.cost_u;
    cost_y = s.cost_y;
    u0 = s.u0;
    y0 = s.y0;
    e0 = s.e0;
    e1 = s.e1;
    u_cl = s.u_cl;
    y_cl = s.y_cl;
    Lf = s.Lf;
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