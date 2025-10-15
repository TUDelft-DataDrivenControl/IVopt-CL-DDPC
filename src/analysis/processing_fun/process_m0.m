function [m0_Uf_mean_iX, m0_Yf_mean_iX, m0_Uf_median_iX, m0_Yf_median_iX, m0_Uf_pctiles_iX, m0_Yf_pctiles_iX] = process_m0(Uf_ivs,Yf_ivs,iX,nu,ny,f,N,seeds,spX,pctiles)
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
% Uf_ivs = {'iv1','iv2a','iv2c','iv3c','iv4a','iv4c','iv5a','iv5c','iv6c'};
% Yf_ivs = {'iv2b','iv3a','iv4b','iv5b','iv6a'};
num_Uf_ivs = numel(Uf_ivs); % needed for nested for loop inside parfor
num_Yf_ivs = numel(Yf_ivs); % needed for nested for loop inside parfor

% initialize cell arrays
m0_Uf_mean_iX = cell(num_Uf_ivs,1);    % mean         cell sizes: nu x ndiags
m0_Yf_mean_iX = cell(num_Yf_ivs,1);    %                          ny x ndiags
m0_Uf_median_iX = m0_Uf_mean_iX;       % median
m0_Yf_median_iX = m0_Yf_mean_iX;
m0_Uf_pctiles_iX = cell(num_Uf_ivs,1); % percentiles  cell sizes: nu x ndiags x num_pctiles
m0_Yf_pctiles_iX = cell(num_Yf_ivs,1); %                          ny x ndiags x num_pctiles

iyf = nu*f + (1:ny*f);
num_diags = f+N-1; % number of anti-diagonals in Uf/Yf (2D case)

%% ---------------- iterating over Uf & Yf -----------------------------
fprintf('\tm0: Uf & Yf IVs done:\t\t');
m0_YUf_tictoc = tic;
for kIVuy = 1:max(num_Uf_ivs,num_Yf_ivs)
    if kIVuy <= num_Uf_ivs
        calc_u = true;
        Uf_iv_name = Uf_ivs{kIVuy};
        m0_Uf_data2 = zeros(nu*f,N,spX); % initialize
    else
        calc_u = false;
    end
    if kIVuy <= num_Yf_ivs
        calc_y = true;
        Yf_iv_name = Yf_ivs{kIVuy};
        m0_Yf_data2 = zeros(ny*f,N,spX); % initialize
    else
        calc_y = false;
    end
    
    % iterate over seeds: get Uf & Yf IV values for all
    % -> parfor is inner loop here to limit data usage (consider size of
    %    Uf_iv for muliple types of IVs)
    parfor ks = 1:spX
        seed = seeds(ks,iX);
        fndata = sprintf('seed_%d.mat',seed);
        Z = load(fndata,'Z').Z;
        
        if calc_u
            Uf_iv = Z.([Uf_iv_name,'_']);
            m0_Uf_data2(:,:,ks) = Uf_iv;
        end
        if calc_y
            switch Yf_iv_name
                case 'iv3a' % IV_Theta: only possible because ny = nlcf (see get_Z.m)
                    Yf_iv = Z.iv3a_(1:ny*f,:);

                case 'iv6a' % Rf_yr0 (future references)
                    Yf_iv = Z.iv6a_;

                otherwise
                    Yf_iv = Z.([Yf_iv_name,'_'])(iyf,:);
            end
            m0_Yf_data2(:,:,ks) = Yf_iv;
        end
    end
    
    %% ---------------------------- for Uf ---------------------------------
    if calc_u
        % get all values belonging to anti-diagonals over seeds
        Uf_ij_adiags = get_subind_diags(f,N,nr=nu,anti=true); % i & j indices for 2D case

        % calculate statistics by iterating over diagonals
        [m0_Uf_mean2, m0_Uf_median2, m0_Uf_pctiles2] = calc_stats(m0_Uf_data2,nu,spX,pctiles,num_diags,Uf_ij_adiags);

        % assign cells in data arrays
        m0_Uf_mean_iX{kIVuy,1}    = m0_Uf_mean2;
        m0_Uf_median_iX{kIVuy,1}  = m0_Uf_median2;
        m0_Uf_pctiles_iX{kIVuy,1} = m0_Uf_pctiles2;
    end

    %% ---------------------------- for Yf ---------------------------------
    if calc_y
        % get all values belonging to anti-diagonals over seeds
        Yf_ij_adiags = get_subind_diags(f,N,nr=ny,anti=true); % i & j indices for 2D case

        % calculate statistics by iterating over diagonals
        [m0_Yf_mean2, m0_Yf_median2, m0_Yf_pctiles2] = calc_stats(m0_Yf_data2,ny,spX,pctiles,num_diags,Yf_ij_adiags);

        % assign cells in data arrays
        m0_Yf_mean_iX{kIVuy,1}    = m0_Yf_mean2;
        m0_Yf_median_iX{kIVuy,1}  = m0_Yf_median2;
        m0_Yf_pctiles_iX{kIVuy,1} = m0_Yf_pctiles2;
    end
end
m0_YUf_time = toc(m0_YUf_tictoc);
fprintf('\t\tFinished in %.2f seconds\n', m0_YUf_time);

end

%% Helper function to calculate statistics for m0
function [m0_Uf_mean2, m0_Uf_median2, m0_Uf_pctiles2] = calc_stats(m0_Uf_data2,nu,spX,pctiles,num_diags,ij_adiags)
% calculates mean, median, percentiles for Uf or Yf data
% written below for Uf data, but used likewise for Yf data

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