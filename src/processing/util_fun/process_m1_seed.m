function m1_UYf_col = process_m1_seed(Uf_ivs, Z, nu, ny, f)
%PROCESS_M1_SEED Compute m1 metrics for a single seed
%   Returns concatenated column vectors m1_Uf_col (numel(Uf_ivs) x 1)
num_Uf_ivs = numel(Uf_ivs);

m1_UYf_col = zeros(num_Uf_ivs,1);

% ---- Uf IVs ----
for kIVu = 1:num_Uf_ivs
    iv_name = Uf_ivs{kIVu};
    switch iv_name
        case 'iv2'
            % leave zero since this is the optimal IV for Uf
        case 'iv3' % IV_Theta: only possible because nu = nlcf (see get_Z.m)
            Uf_iv = Z.iv3_(1:ny*f,:);
            m1_UYf_col(kIVu) = norm(Uf_iv - Z.iv2_, 'fro');
        otherwise
            Uf_iv = Z.([iv_name,'_']);
            m1_UYf_col(kIVu) = norm(Uf_iv - Z.iv2_, 'fro');
    end
end

end
