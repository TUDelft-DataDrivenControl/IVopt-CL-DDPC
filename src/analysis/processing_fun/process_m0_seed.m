function UYf_row = process_m0_seed(Z, Uf_ivs, Yf_ivs, nu, ny, f)
%PROCESS_M0_SEED Load one seed file and extract Uf and Yf IV matrices
%   UYf_row is a 1 x (num_Uf_ivs+num_Yf_ivs) cell array where the first
%   num_Uf_ivs cells are Uf matrices and the remaining are Yf matrices.

    num_Uf_ivs = numel(Uf_ivs);
    num_Yf_ivs = numel(Yf_ivs);

    Uf_temp = cell(1, num_Uf_ivs);
    Yf_temp = cell(1, num_Yf_ivs);

    iyf = nu*f + (1:ny*f);

    % Uf IVs
    for kIVuy = 1:num_Uf_ivs
        Uf_iv_name = Uf_ivs{kIVuy};
        Uf_iv = Z.([Uf_iv_name,'_']);
        Uf_temp{kIVuy} = Uf_iv;
    end

    % Yf IVs
    for kIVuy = 1:num_Yf_ivs
        Yf_iv_name = Yf_ivs{kIVuy};
        switch Yf_iv_name
            case 'iv3a'
                Yf_iv = Z.iv3a_(1:ny*f,:); % only select \Xi_f (IV_Theta). N.B. works here because size(\Xi_f,1)/f==ny
            case 'iv6a'
                Yf_iv = Z.iv6a_;
            otherwise
                Yf_iv = Z.([Yf_iv_name,'_'])(iyf,:);
        end
        Yf_temp{kIVuy} = Yf_iv;
    end

    UYf_row = [Uf_temp Yf_temp];
end
