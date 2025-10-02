function [m1_Uf_data_iX, m1_Yf_data_iX, m2_data_iX, m3_data_iX] = process_m123(Uf_ivs,Yf_ivs,Cases,IDerrorTypes,cost_types,iX,nu,ny,f,seeds,spX)

% initializing max counters needed for nested for loop inside parfor
num_Uf_ivs       = numel(Uf_ivs); 
num_Yf_ivs       = numel(Yf_ivs); 
num_Cases        = numel(Cases);
num_IDerrorTypes = numel(IDerrorTypes);
num_cost_types   = numel(cost_types);

m1_Uf_data_iX = zeros(num_Uf_ivs,spX);                    % cell sizes: num_Uf_ivs x spX
m1_Yf_data_iX = zeros(num_Yf_ivs,spX);                    % cell sizes: num_Yf_ivs x spX
m2_data_iX    = zeros(num_IDerrorTypes,num_Cases,spX);    % cell sizes: num_IDerrorTypes x num_Cases x spX
m3_data_iX    = zeros(num_cost_types,num_Cases,spX);      % cell sizes: num_cost_types   x num_Cases x spX

iyf = nu*f + (1:ny*f); % row indices of Yf in Z

% Ensure we have a parallel pool
if isempty(gcp('nocreate'))
    parpool;
end

% Initialize progress tracking
w = waitbar(0,sprintf('Progress: (%d/%d) -> %.1f%%',0,spX,0));
D = parallel.pool.DataQueue;
afterEach(D,@parforWaitbar);
parforWaitbar(w, spX); % Initialize waitbar function

m123_tictoc = tic;
txt_iters = sprintf('\tm1, m2, m3: Iterating over seeds (see waitbar).');
fprintf('%s',txt_iters);
parfor ks = 1:spX
    seed = seeds(ks,iX);
    fndata = sprintf('seed_%d.mat',seed);
    [FroIDerror,Z,cost_tot,cost_u,cost_y] = load_seedmat(fndata);
    
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
    
    send(D, []);
end % of parfor
m123_time = toc(m123_tictoc);
close(w);
fprintf([repmat('\b',1,numel(txt_iters)),'\tm1, m2, m3: Iterating over seeds\t\tFinished in %.2f seconds\n'], m123_time);

end

%% Helper functions
function [FroIDerror,Z,cost_tot,cost_u,cost_y] = load_seedmat(fnpath)
    % Only load necessary variables to improve memory usage and performance
    s = load(fnpath,'FroIDerror','Z','cost_tot','cost_u','cost_y');
    FroIDerror = s.FroIDerror;
    Z = s.Z;
    cost_tot = s.cost_tot;
    cost_u = s.cost_u;
    cost_y = s.cost_y;
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