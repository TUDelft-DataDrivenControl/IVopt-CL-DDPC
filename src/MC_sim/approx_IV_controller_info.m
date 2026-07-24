function [Uf_iv5a, Uf_iv5b, Uf_iv5c, Uf_iv5d] = approx_IV_controller_info(u,y,w,opts,Cz)
%% Approximates the optimal IVs Uf using controller information
% This function approximates the IV matrix Uf that would have
% been obtained without future noise Ef, using knowledge of the
% functional form of the 1 DOF feedback controller C_{fb}:
% 
% feedback law:
%       u_k = C_{fb}(q) (w_k - y_k)
%
% Resulting composition of data matrix D:
%       D = [U_{\varrho}; Y_{\varrho}; W_{\nu}; W_f]

[gamma,p,f,DMCS] = deal(opts.gamma,opts.p,opts.f,opts.DMCS);
[nu,Nbar] = size(u);
ny = size(y,1);
validateattributes(y, {'double'},{'size',[ny Nbar]});
validateattributes(w, {'double'},{'size',[ny Nbar]});

%% determine past data length to use - varrho
varrho = max(gamma,p); % get max of past controller & plant windows

%% ===== make the actual matrices pertaining to the controller ============

[Ac,Bc,Cc,Dc] = ssdata(Cz);
[Lce, Lcu, Tcf] = get_lifted_mats(Ac,Bc,Cc,Dc,[],gamma,f);

%% ===================== estimating system matrices =======================

% make data matrices matrices
Uvf = make_Page(u,varrho+f,DMCS); Uv = Uvf(1:nu*varrho,:); Uf = Uvf(nu*varrho+1:end,:);   
Yvf = make_Page(y,varrho+f,DMCS); Yv = Yvf(1:ny*varrho,:); Yf = Yvf(ny*varrho+1:end,:);
Wvf = make_Page(w,varrho+f,DMCS); Wr = Wvf(1:ny*varrho,:); Wf = Wvf(ny*varrho+1:end,:); % Wr -> gamma instead of varrho deep

% select specific gamma and p long matrices
Wr = Wr(end-gamma*ny+1:end,:);
Up = Uv(end-p*nu+1:end,:); Ur = Uv(end-gamma*nu+1:end,:);
Yp = Yv(end-p*ny+1:end,:); Yr = Yv(end-gamma*ny+1:end,:);

% contribibution from Wr and Wf to rf
Rw_rf = [Lce Tcf]*[Wr;Wf] + Lcu*Ur - Lce*Yr; % f*nu x N

DMr = [Up;Yp;Rw_rf];


iv_names = {'iv5a','iv5b','iv5c','iv5d'};
iv_names = intersect(iv_names,opts.Cases);
[Uf_iv5a,Uf_iv5b,Uf_iv5c,Uf_iv5d] = deal(nan); % predefine for output
for kiv = 1:length(iv_names)
iv_name = iv_names{kiv};
switch iv_name
%% ------------------------------ IV5a (2SLS) -----------------------------
% NOTE: this is the algorithm reported as `IV5' in the paper
case 'iv5a'
L_clu = Uf*pinv(DMr); % -> estimates [Guu Guy Mu; Gyu Gyy Tuf*Mu]
Uf_iv5a = L_clu*DMr;

%% ------------------------------ IV5b (2SLS + causal + time invariant) ---
case 'iv5b'
% -> attempt to improve L_clu using causality & time invariance
L_clu1 = L_clu;

% part for Wf should be causal -> impose block-tril structure
L_clu1(:,end-f*nu+1:end) = make_blk_tril(L_clu1(:,end-f*nu+1:end),[nu,nu]);

% part for Wf should be time invariant -> block-diagonal averaging
L_clu1(:,end-f*nu+1:end) = blk_toeplitz_mean(L_clu1(:,end-f*nu+1:end),nu,nu);

Uf_iv5b = L_clu1*DMr;

%% ----------------------------- IV5c (row-by-row) -------------------------
case 'iv5c'
% -> estimate L_clu row by row (-> causal by construction)
L_clu2 = zeros(nu*f,size(DMr,1));
rows = 1:nu;
colend = (nu+ny)*p + nu;
for k = 1:f
    L_clu2(rows,1:colend) = Uf(rows,:)*pinv(DMr(1:colend,:));
    rows = rows + nu;
    colend = colend + nu;
end

Uf_iv5c = L_clu2*DMr;

%% ----------------------------- IV4d (row-by-row + time invariant) -----------
case 'iv5d'
L_clu3 = L_clu2; % causal by construction
% part for Wf should be time invariant -> block-diagonal averaging
L_clu3(:,end-f*nu+1:end) = blk_toeplitz_mean(L_clu3(:,end-f*nu+1:end),nu,nu);

Uf_iv5d = L_clu3*DMr;

end % of switch
end % of for loop
end

%% Local functions
function [Lup,Lyp,Tuf, varargout] = get_lifted_mats(A,B,C,D,K,p,f)
    nyC = size(C,1);
    if isempty(K)
        nxA = size(A,1);
        K = zeros(nxA,nyC);
    end

    % make modified controllability matrices
    tA = A-K*C;
    tB = B-K*D;
    tAp_pinv = tA^p*pinv(make_ext_obsv(A,C,p));
    Kyp = make_ext_ctrb(tA,K, p,rev=true) + tAp_pinv;
    Kup = make_ext_ctrb(tA,tB,p,rev=true) - tAp_pinv*make_blk_tril_toeplitz(A,B,C,D,p);

    Gamf = make_ext_obsv(A,C,f);

    Lup = Gamf*Kup;
    Lyp = Gamf*Kyp;
    Tuf = make_blk_tril_toeplitz(A,B,C,D,f);
    
    if nargout >= 4
        varargout{1}  = make_blk_tril_toeplitz(A,K,C,eye(nyC),f);
    end
    if nargout >= 5
        Hp  = make_blk_tril_toeplitz(A,K,C,eye(nyC),p);
        varargout{2}  = Gamf*tAp_pinv*Hp;
    end
end

