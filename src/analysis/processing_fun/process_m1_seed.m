function m1_UYf_col = process_m1_seed(Uf_ivs, Yf_ivs, Z, nu, ny, f)
%PROCESS_M1_SEED Compute m1 metrics for a single seed
%   Returns concatenated column vectors m1_Uf_col (numel(Uf_ivs) x 1) and
%   m1_Yf_col (numel(Yf_ivs) x 1) for the supplied Z struct.

num_Uf_ivs = numel(Uf_ivs);
num_Yf_ivs = numel(Yf_ivs);

m1_UYf_col = zeros(num_Uf_ivs+num_Yf_ivs,1);

% precompute index for Yf rows inside Z when needed
iyf = nu*f + (1:ny*f);

% ---- Uf IVs ----
for kIVu = 1:num_Uf_ivs
    iv_name = Uf_ivs{kIVu};
    switch iv_name
        case 'iv2a'
            % leave zero since this is the optimal IV for Uf
        otherwise
            Uf_iv = Z.([iv_name,'_']);
            m1_UYf_col(kIVu) = norm(Uf_iv - Z.iv2a_, 'fro');
    end
end

% ---- Yf IVs ----
for kIVy = 1:num_Yf_ivs
    iv_name = Yf_ivs{kIVy};

    switch iv_name
        case 'iv2b'
            % leave zero since this contains the optimal IV for Yf
            % set calc_norm = false to keep zero and suppress warnings
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
        m1_UYf_col(kIVy + num_Uf_ivs) = norm(Yf_iv - Z.iv2b_(iyf,:), 'fro');
    end
end
end
