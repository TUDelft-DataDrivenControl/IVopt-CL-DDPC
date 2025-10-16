function [m0_Uf_mean_iX, m0_Yf_mean_iX, m0_Uf_median_iX, m0_Yf_median_iX, m0_Uf_pctiles_iX, m0_Yf_pctiles_iX] = process_m0(Uf_ivs,Yf_ivs,iX,nu,ny,f,N,seeds,spX,pctiles)
% ======================== initialize measure 0 (m0) ======================
% Restructured version with parfor over seeds (ks) and inner loop over IVs

num_Uf_ivs = numel(Uf_ivs);
num_Yf_ivs = numel(Yf_ivs);
num_IVs = num_Uf_ivs + num_Yf_ivs;

% initialize cell arrays for final results
m0_UYf_mean_iX = cell(num_IVs,1);
m0_UYf_median_iX = m0_UYf_mean_iX;
m0_UYf_pctiles_iX = cell(num_IVs,1);

iyf = nu*f + (1:ny*f);
num_diags = f+N-1;

% Pre-compute diagonal indices (same for all seeds)
Uf_ij_adiags = get_subind_diags(f,N,nr=nu,anti=true);
Yf_ij_adiags = get_subind_diags(f,N,nr=ny,anti=true);

%% ---------------- Parallel iteration over seeds -------------------------
fprintf('\tm0: Processing seeds in parallel...\t\t');
m0_YUf_tictoc = tic;

% Initialize cell arrays to collect data from each seed
% Structure: {ks}{kIVuy} = data for seed ks, IV kIVuy
UYf_data_per_seed = cell(spX, num_IVs);

parfor ks = 1:spX
    seed = seeds(ks,iX);
    fndata = sprintf('seed_%d.mat',seed);
    Z = load(fndata,'Z').Z;
    
    % Temporary storage for this seed
    Uf_temp = cell(1,num_Uf_ivs);
    Yf_temp = cell(1,num_Yf_ivs);
    
    % Inner loop over Uf IVs
    for kIVuy = 1:num_Uf_ivs
        Uf_iv_name = Uf_ivs{kIVuy};
        Uf_iv = Z.([Uf_iv_name,'_']);
        Uf_temp{kIVuy} = Uf_iv;  % Store (nu*f x N) matrix
    end
    
    % Inner loop over Yf IVs
    for kIVuy = 1:num_Yf_ivs
        Yf_iv_name = Yf_ivs{kIVuy};
        
        switch Yf_iv_name
            case 'iv3a' % IV_Theta
                Yf_iv = Z.iv3a_(1:ny*f,:);
            case 'iv6a' % Rf_yr0
                Yf_iv = Z.iv6a_;
            otherwise
                Yf_iv = Z.([Yf_iv_name,'_'])(iyf,:);
        end
        Yf_temp{kIVuy} = Yf_iv;  % Store (ny*f x N) matrix
    end
    
    % Assign to collection arrays
    UYf_data_per_seed(ks, :) = [Uf_temp Yf_temp];
end

%% ---------------- Aggregate and calculate statistics --------------------
fprintf('Done.\n\tm0: Calculating statistics...\t\t');

parfor kIVuy = 1:num_IVs
    if kIVuy <= num_Uf_ivs
        % ==== Process Uf ====
        m0_data3 = zeros(nu*f, N, spX);
        Uf_data_per_seed_kIV = UYf_data_per_seed(:,kIVuy)
        for ks = 1:spX
            m0_data3(:,:,ks) = Uf_data_per_seed_kIV{ks, 1};
        end
        
        [m0_mean2, m0_median2, m0_pctiles2] = ...
            calc_stats(m0_data3, nu, spX, pctiles, num_diags, Uf_ij_adiags);
        
    else
        % ==== Process Yf ====
        % kY = kIVuy - num_Uf_ivs; % adjust index for Yf section

        m0_data3 = zeros(ny*f, N, spX);
        Yf_data_per_seed_kY = UYf_data_per_seed(:,kIVuy);
        for ks = 1:spX
            m0_data3(:,:,ks) = Yf_data_per_seed_kY{ks, 1};
        end
        
        [m0_mean2, m0_median2, m0_pctiles2] = ...
            calc_stats(m0_data3, ny, spX, pctiles, num_diags, Yf_ij_adiags);
        
    end
    m0_UYf_mean_iX{kIVuy,1}    = m0_mean2;
    m0_UYf_median_iX{kIVuy,1}  = m0_median2;
    m0_UYf_pctiles_iX{kIVuy,1} = m0_pctiles2;
end
m0_Uf_mean_iX    = m0_UYf_mean_iX(1:num_Uf_ivs,1);
m0_Yf_mean_iX    = m0_UYf_mean_iX(num_Uf_ivs+1:end,1);    clear m0_UYf_mean_iX
m0_Uf_median_iX  = m0_UYf_median_iX(1:num_Uf_ivs,1);
m0_Yf_median_iX  = m0_UYf_median_iX(num_Uf_ivs+1:end,1);  clear m0_UYf_median_iX
m0_Uf_pctiles_iX = m0_UYf_pctiles_iX(1:num_Uf_ivs,1);
m0_Yf_pctiles_iX = m0_UYf_pctiles_iX(num_Uf_ivs+1:end,1); clear m0_UYf_pctiles_iX

m0_YUf_time = toc(m0_YUf_tictoc);
fprintf('Finished in %.2f seconds\n', m0_YUf_time);

end

%% Helper function to calculate statistics for m0
function [m0_Uf_mean2, m0_Uf_median2, m0_Uf_pctiles2] = calc_stats(m0_Uf_data2,nu,spX,pctiles,num_diags,ij_adiags)
% calculates mean, median, percentiles for Uf or Yf data

num_pctiles = numel(pctiles);
% initialize m0_Uf_...
m0_Uf_mean2    = zeros(nu,num_diags);             % mean
m0_Uf_median2  = m0_Uf_mean2;                     % median
m0_Uf_pctiles2 = zeros(nu,num_diags,num_pctiles); % percentiles

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

end