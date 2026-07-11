function [Uf_iv4a,Uf_iv4b,Uf_iv4c,Uf_iv4d] = approx_IV_no_controller_info(u,y,w,opts)
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

iv_names = {'iv4a','iv4b','iv4c','iv4d'};
iv_names = intersect(iv_names,opts.Cases);
[Uf_iv4a,Uf_iv4b,Uf_iv4c,Uf_iv4d] = deal(nan); % predefine for output
for kiv = 1:length(iv_names)
iv_name = iv_names{kiv};
switch iv_name
%% ----------------------------- (2SLS) -----------------------------------
case 'iv4a' 
L_clu = Uf*pinv(D); % initial estimate
Uf_iv4a = L_clu*D;

%% ----------------------- (2SLS + causal + time invariant) ---------------
case 'iv4b'
% -> attempt to improve L_clu using causality & time invariance
L_clu1 = L_clu;

% part for Wf should be causal -> impose block-tril structure
L_clu1(:,end-f*ny+1:end) = make_blk_tril(L_clu1(:,end-f*ny+1:end),[nu,ny]);

% part for Wf should be time invariant -> block-diagonal averaging
L_clu1(:,end-f*ny+1:end) = blk_toeplitz_mean(L_clu1(:,end-f*ny+1:end),nu,ny);

Uf_iv4b = L_clu1*D;

%% ------------------------- (row-by-row) ---------------------------------
case 'iv4c'
% -> estimate L_clu row by row (-> causal by construction)
L_clu2 = zeros(nu*f,size(D,1));
rows = 1:nu;
colend = (nu+ny)*varrho + ny*gamma + ny;
for k = 1:f
    L_clu2(rows,1:colend) = Uf(rows,:)*pinv(D(1:colend,:));
    rows = rows + nu;
    colend = colend + ny;
end

Uf_iv4c = L_clu2*D;

%% ---------------------- (row-by-row + time invariant) -------------------
case 'iv4d'
L_clu3 = L_clu2; % causal by construction
% part for Wf should be time invariant -> block-diagonal averaging
L_clu3(:,end-f*ny+1:end) = blk_toeplitz_mean(L_clu3(:,end-f*ny+1:end),nu,ny);

Uf_iv4d = L_clu3*D;

end % of switch
end % of for loop

end