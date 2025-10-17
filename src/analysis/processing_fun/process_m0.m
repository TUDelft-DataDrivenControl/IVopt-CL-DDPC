function m0_UYf_iX = process_m0(Uf_ivs,Yf_ivs,iX,nu,ny,f,N,seeds,spX,pctiles)
% ======================== initialize measure 0 (m0) ======================
% Restructured version with parfor over seeds (ks) and inner loop over IVs

num_Uf_ivs = numel(Uf_ivs);
num_Yf_ivs = numel(Yf_ivs);
num_IVs = num_Uf_ivs + num_Yf_ivs;

% results will be produced after aggregation

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
    Z = load(fndata, 'Z').Z;
    UYf_data_per_seed(ks, :) = process_m0_seed(Z, Uf_ivs, Yf_ivs, nu, ny, f);
end

%% ---------------- Aggregate and calculate statistics --------------------
fprintf('Done.\n\tm0: Calculating statistics...\t\t');

m0_UYf_iX = calc_stats_seeds_m0(UYf_data_per_seed, num_Uf_ivs, nu, ny, f, N, spX, pctiles, Uf_ij_adiags, Yf_ij_adiags);

m0_YUf_time = toc(m0_YUf_tictoc);
fprintf('Finished in %.2f seconds\n', m0_YUf_time);

end

% aggregation moved to calc_stats_seeds_m0.m