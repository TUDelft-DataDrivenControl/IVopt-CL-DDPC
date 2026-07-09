function [Uf_2sls,Uf_iv4] = approx_IV_no_controller_info(u,y,w,opts)
%% Approximates the optimal IVs Uf & Yf without controller information (Algorithm 1 from the article)
% This function approximates the IV matrices Uf and Yf that would have
% been obtained without future noise Ef, using no knowledge of the
% functional form of the 1 DOF feedback controller C_{fb}:
% 
% feedback law:
%       u_k = C_{fb}(q) (w_k - y_k)
%
% The only knowledge that is assumed of C_{fb} is an accurate upper bound
% of its lag: \gamma
%
% Resulting composition of data matrix D:
%       D = [U_{\varrho}; Y_{\varrho}; W_{\nu}; W_f]

[gamma,p,f,nu,ny,DMCS] = deal(opts.gamma,opts.p,opts.f,opts.nu,opts.ny,opts.DMCS);
varrho = max(gamma,p);

[~,Nbar] = size(u);
validateattributes(u, {'double'},{'size',[nu Nbar]});
validateattributes(y, {'double'},{'size',[ny Nbar]});
validateattributes(w, {'double'},{'size',[ny Nbar]});

% create data matrix D and Uf
Uvf = make_Page(u,varrho+f,DMCS); Uv = Uvf(1:nu*varrho,:); Uf = Uvf(nu*varrho+1:end,:);   
Yvf = make_Page(y,varrho+f,DMCS); Yv = Yvf(1:ny*varrho,:);
Wvf = make_Page(w,varrho+f,DMCS); Wr = Wvf(1:ny*varrho,:); Wf = Wvf(ny*varrho+1:end,:);
Wr = Wr(end-gamma+1:end,:); % select last gamma block rows
D = [Uv;Yv;Wr;Wf];

%% ------------------------------ 2SLS (IV7) ------------------------------
L_clu = Uf*pinv(D); % initial estimate
Uf_2sls = L_clu*D;

%% ------------------------------ IV4a ------------------------------------
% -> attempt to improve L_clu using causality & time invariance
L_clu1 = L_clu;
% part for Wf should be causal & time invariant -> impose lower (block-)tril (averaged) structure
L_clu1(:,end-f*ny+1:end) = blk_tril_avg(L_clu1(:,end-f*ny+1:end),nu,ny);

Uf_iv4 = L_clu1*D;

end

function Lpart = blk_tril_avg(Lpart,nrow,ncol)
[s1,s2] = size(Lpart);
blkrows = s1/nrow;
blkcols = s2/ncol; 
if blkrows~=floor(blkrows) || blkcols~=floor(blkcols)
    error('The block dimensions do not fit evenly in the provided matrix.');
elseif blkrows ~= blkcols
    error('The number of block rows must be equal to the number of block columns');
end
% part for Wf should be causal -> impose (block-)tril structure
Lpart = Lpart.*kron(tril(ones(blkrows)),ones(nrow,ncol)); % blocks: nrow x ncol

% take mean of corresponding elements in blocks along diagonals
Lpart = blk_toeplitz_mean(Lpart,nrow,ncol);
end