function [Uf_iv,Yf_iv] = approx_IV_controller_info(u,y,w,opts,Cz)
%% Approximates the optimal IVs Uf & Yf using controller information (Algorithm 2 from the article)
% This function approximates the IV matrices Uf and Yf that would have
% been obtained without future noise Ef, using knowledge of the
% functional form of the 1 DOF feedback controller C_{fb}:
% 
% feedback law:
%       u_k = C_{fb}(q) (w_k - y_k)
%
% Resulting composition of data matrix D:
%       D = [U_{\varrho}; Y_{\varrho}; W_{\nu}; W_f]

[gamma,p,f] = deal(opts.gamma,opts.p,opts.f);
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

% make hankel matrices
[~,Uv,Uf] = make_Hankel(u,varrho,f);   
[~,Yv,Yf] = make_Hankel(y,varrho,f); 
[~,Wr,Wf] = make_Hankel(w,varrho,f); % Wr -> gamma instead of varrho deep

% select specific gamma and p long matrices
Wr = Wr(end-gamma*ny+1:end,:);
Up = Uv(end-p*nu+1:end,:); Ur = Uv(end-gamma*nu+1:end,:);
Yp = Yv(end-p*ny+1:end,:); Yr = Yv(end-gamma*ny+1:end,:);

% contribibution from Wr and Wf to rf
Rw_rf = [Lce Tcf]*[Wr;Wf] + Lcu*Ur - Lce*Yr; % f*nu x N

DMr = [Up;Yp;Rw_rf];

iuf = 1:nu*f;
iyf = f*nu+(1:ny*f);
pnuy = p*(nu+ny);

%% initial estimate:
G_0 = [Uf;Yf]*pinv(DMr); % -> estimates [Guu Guy Mu; Gyu Gyy Tuf*Mu]

% extract Mu & Tuf*Mu estimates
Mu_0    = G_0(iuf,pnuy+(1:f*nu));   % Mu
TufMu_0 = G_0(iyf,pnuy+(1:f*nu));   % Tuf*Mu

% block-lower triangularize & average
Mu_1    = blk_tril_avg(Mu_0,   nu,nu);
TufMu_1 = blk_tril_avg(TufMu_0,ny,nu);

% re-estimate remaining parameters taking into account contribution of [Mu_1;TufMu_1]*Rw_rf
G_1 = ([Uf;Yf]-[Mu_1;TufMu_1]*Rw_rf)*pinv([Up;Yp]);

%% estimate Uf_iv & Yf_iv
UY_iv = [G_1 [Mu_1;TufMu_1]]*DMr;
Uf_iv = UY_iv(iuf,:);
Yf_iv = UY_iv(iyf,:);

end

%% Local functions
function [Lup,Lyp,Tuf, varargout] = get_lifted_mats(A,B,C,D,K,p,f)
    nyC = size(C,1);
    if isempty(K)
        nxA = size(A,1);
        K = zeros(nxA,nyC);
    end

    % make 'effective' controllability matrices
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

function Mat = blk_tril_avg(Mat,dim1,dim2)
    Mat = make_blk_tril(Mat, [dim1,dim2]);
    Mat = blk_toeplitz_mean(Mat,dim1,dim2);
end

