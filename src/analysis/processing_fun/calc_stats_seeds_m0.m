function m0_UYf_iX = calc_stats_seeds_m0(UYf_data_per_seed, num_Uf_ivs, nu, ny, f, N, spX, pctiles)
%CALC_STATS_SEEDS_M0 Aggregate per-seed U/Yf data and compute stats per IV
%   Inputs:
%     UYf_data_per_seed : cell(spX, num_IVs) with each cell a matrix (nu*f x N) or (ny*f x N)
%     num_Uf_ivs        : number of Uf IVs (first block of columns)
%     nu, ny, f, N, spX : dimensions
%     pctiles           : percentiles vector
%     Uf_ij_adiags, Yf_ij_adiags : cell arrays with diag indices

% Pre-compute cell arrays diagonal indices (same for all seeds)
Uf_ij_adiags = get_subind_diags(f,N,nr=nu,anti=true);
Yf_ij_adiags = get_subind_diags(f,N,nr=ny,anti=true);

num_IVs = size(UYf_data_per_seed, 2);
num_diags = f + N - 1;

m0_UYf_iX = cell(num_IVs,3); % (:,1)=mean, (:,2)=median, (:,3)=pctiles

for kIVuy = 1:num_IVs
    if kIVuy <= num_Uf_ivs
        % ==== Process Uf ====
        m0_data3 = zeros(nu*f, N, spX);
        for ks = 1:spX
            m0_data3(:,:,ks) = UYf_data_per_seed{ks, kIVuy};
        end
        [m0_UYf_iX{kIVuy,1}, m0_UYf_iX{kIVuy,2}, m0_UYf_iX{kIVuy,3}] = calc_stats_inner(m0_data3, nu, spX, pctiles, num_diags, Uf_ij_adiags);
    else
        % ==== Process Yf ====
        m0_data3 = zeros(ny*f, N, spX);
        for ks = 1:spX
            m0_data3(:,:,ks) = UYf_data_per_seed{ks, kIVuy};
        end
        [m0_UYf_iX{kIVuy,1}, m0_UYf_iX{kIVuy,2}, m0_UYf_iX{kIVuy,3}] = calc_stats_inner(m0_data3, ny, spX, pctiles, num_diags, Yf_ij_adiags);
    end
end

end

function [m0_mean2, m0_median2, m0_pctiles2] = calc_stats_inner(m0_data3, n, spX, pctiles, num_diags, ij_adiags)
%CALC_STATS_INNER Compute mean/median/pctiles for stacked 3D data
num_pctiles = numel(pctiles);
m0_mean2    = zeros(n, num_diags);
m0_median2  = m0_mean2;
m0_pctiles2 = zeros(n, num_diags, num_pctiles);

for kd = 1:num_diags
    rows = ij_adiags{kd}(:,1);
    cols = ij_adiags{kd}(:,2);
    rows2 = repmat(rows, spX, 1);
    cols2 = repmat(cols, spX, 1);
    i3s = kron((1:spX).', ones(numel(cols),1));
    idxlin = sub2ind(size(m0_data3), rows2, cols2, i3s);
    Uf_sel = m0_data3(idxlin);
    Uf_sel = reshape(Uf_sel, n, []);
    m0_mean2(:,kd)      = mean(Uf_sel,2);
    m0_median2(:,kd)    = median(Uf_sel,2);
    m0_pctiles2(:,kd,:) = prctile(Uf_sel, pctiles, 2);
end
end
