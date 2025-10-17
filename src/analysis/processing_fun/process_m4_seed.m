function [Yf_RelErr_sd_block, Yf_RelErr_mean_block] = process_m4_seed(e0,e1,u0,y0,u_cl,y_cl,Lf,effEpMat,Hf,Cases,p,f)
%PROCESS_M4_SEED Compute Yf relative error statistics for a single seed
%   Returns two blocks of size (ny*f) x num_Cases x 4 where the 4 corresponds
%   to the different error components. The function assumes sizes are
%   consistent with caller; it computes per-case statistics and returns
%   matrices that the caller can place into the per-seed arrays.

% build noise contributions
e_all = [e0(:,end-p+1:end) e1];
[~, Ep, Ef] = make_Hankel(e_all, p, f);
Yf_by_Ep = effEpMat*Ep; % past noise contribution to Yf
Yf_by_Ef = Hf*Ef;       % future noise contribution to Yf

num_Cases = numel(Cases);

% Determine ny*f from dimensions of Yf_by_Ep
nyf = size(Yf_by_Ep,1);

Yf_RelErr_block = zeros(nyf, num_Cases, 4, 2); % (:,:,:,1) -> std. dev, (:,:,:,2) -> mean

for kC = 1:num_Cases
    caseName = Cases{kC};
    % create Hankel matrices for inputs & outputs
    u_all = [u0(:,end-p+1:end) u_cl.(caseName)];
    [~, Up, Uf] = make_Hankel(u_all, p, f);
    y_all = [y0(:,end-p+1:end) y_cl.(caseName)];
    [~, Yp, Yf] = make_Hankel(y_all, p, f);

    Yf_hat = Lf.(caseName)*[Up;Yp;Uf];      % prediction w/ Lf estimate
    if ~strcmp(caseName,'actLf')
        Yf_hatS = Lf.('actLf')*[Up;Yp;Uf];  % prediction w/ actual Lf
    else
        Yf_hatS = Yf_hat;
    end

    % since data can vary greatly in magnitude: use relevative difference
    [Yf_RelErr_block(:,kC,1,1), Yf_RelErr_block(:,kC,1,2)] = std( (Yf_hat  - Yf)./Yf, 0, 2);
    [Yf_RelErr_block(:,kC,2,1), Yf_RelErr_block(:,kC,2,2)] = std( (Yf_hatS - Yf)./Yf, 0, 2);
    [Yf_RelErr_block(:,kC,3,1), Yf_RelErr_block(:,kC,3,2)] = std( Yf_by_Ep./Yf, 0, 2);
    [Yf_RelErr_block(:,kC,4,1), Yf_RelErr_block(:,kC,4,2)] = std( Yf_by_Ef./Yf, 0, 2);
end

end
