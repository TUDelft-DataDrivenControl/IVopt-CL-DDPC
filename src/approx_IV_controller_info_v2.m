function [Uf_iv,Yf_iv,xCz_iv] = approx_IV_controller_info_v2(u,y,wr,rho,p,f,Ziv_init,Cz,nit_max,plant,e)
% This function approximates the IV matrices Uf and Yf that would have
% been obtained without future noise Ef, using knowledge of the
% functional form of the 1 DOF feedback controller C_{fb}:
% 
% feedback law:
%       u_k = C_{fb}(q) (w_k - y_k)
%
% The only knowledge that is assumed of C_{fb} is an accurate upper bound
% of its lag: rho
%
% Resulting composition of data matrix D:
%       D = [U_{\varrho}; Y_{\varrho}; W_{\rho}; W_f]

[nu,Nbar] = size(u);
ny = size(y,1);
validateattributes(y,  {'double'},{'size',[ny Nbar]});
validateattributes(wr,{'double'},{'size',[ny Nbar]});

%% determine past data length to use - varrho

% nxc = size(Ac,1);
% Gamma_c = make_ext_obsv(Ac,Cc,nxc);
% if rank(Gamma_c)~= nxc
%     error('Controller is not observable');
% end
% % rho = lag of the controller
% rho = nan;
% for k = 1:nxc
%     if rank(Gamma_c(1:nu*k,:)) == nxc
%         rho = k;
%         break;
%     end
% end

% get max of past controller & plant windows
varrho = max(rho,p);

%% ===================== make the actual matrices =========================

function [Ksu,Ksy] = get_Ksuy(A,B,C,D,M,s)
    if isempty(M)
        nA = size(A,1);
        nC = size(C,1);
        M = zeros(nA,nC);
    end
    tA = A-M*C;
    tB = B-M*D;
    tAs_pinv = tA^s*pinv(make_ext_obsv(A,C,s));
    Ksy = make_ext_ctrb(tA,M, s,rev=true) + tAs_pinv;
    Ksu = make_ext_ctrb(tA,tB,s,rev=true) - tAs_pinv*make_blk_tril_toeplitz(A,B,C,D,s);
end

%==========================================================================
%--------------------------- the controller -------------------------------
%==========================================================================
% (1)       u_f = r_f - Lcur u_rho - Lcyr y_rho - Tcf y_f
% where:
%   rho   = lag of the controller
%   u_f   = future inputs             (nu*f   | 1)
%   r_f   = future exogenous inputs   (nu*f   | 1)
%   u_rho = past inputs               (nu*rho | 1)
%   y_rho = past outputs              (ny*rho | 1)
%   y_f   = future outputs            (ny*f   | 1)
[Ac,Bc,Cc,Dc] = ssdata(Cz);
% (extended) observability matrices
Gam_cf = make_ext_obsv(Ac,Cc,f);
Gam_cr = make_ext_obsv(Ac,Cc,rho);

% (extended) reversed controllability matrices
Kcr = make_ext_ctrb(Ac,Bc,rho,rev=true);

% Toeplitz matrices
Tcf = make_blk_tril_toeplitz(Ac,Bc,Cc,Dc,f);
Tcr = make_blk_tril_toeplitz(Ac,Bc,Cc,Dc,rho);

% other matrices
Acr = Ac^rho;
Acrpinv = Acr*pinv(Gam_cr);

% construct Lcur & Lcyr
Lcur = Gam_cf*Acrpinv;
Lcyr = Gam_cf*(Kcr-Acrpinv*Tcr);

% [Lcur,Lcyr] = get_Ksuy(Ac,Bc,Cc,Dc,[],rho);
% Lcur = Gam_cf*Lcur;
% Lcyr = Gam_cf*Lcyr;

%==========================================================================
%------------------------------- the plant --------------------------------
%==========================================================================
% (2)       y_f = Lup u_p + Lyp y_p + Tuf u_f + Hf e_f - Be e_p
% where
%   u_p = past inputs               (nu*p | 1)
%   y_p = past outputs              (ny*p | 1)
%   e_f = future noise              (ny*f | 1)
%   e_p = past noise                (ny*p | 1)

% A, B, C, D, K of actual plant
[A,BK,C,DI] = ssdata(plant);
B = BK(:,1:nu); K = BK(:,end-ny+1:end);
D = DI(:,1:nu);

% (extended) observability matrices
Gam_f = make_ext_obsv(A,C,f);
Gam_p = make_ext_obsv(A,C,p);

% (extended) reversed controllability matrices
tKpu = make_ext_ctrb(A-K*C,B-K*D,p,rev=true);
tKpy = make_ext_ctrb(A-K*C,K,p,rev=true);

% Toeplitz matrices
Hf  = make_blk_tril_toeplitz(A,K,C,eye(ny),f);
Hp  = make_blk_tril_toeplitz(A,K,C,eye(ny),p);
Tuf = make_blk_tril_toeplitz(A,B,C,D,f);            % <-- need to estimate
Tup = make_blk_tril_toeplitz(A,B,C,D,p);

% other matrices
tAp = (A-K*C)^p;
tAppinv = tAp*pinv(Gam_p);

% construct Lup, Lyp, Be
Lup = Gam_f*(tKpu - tAppinv*Tup);
Lyp = Gam_f*(tKpy + tAppinv);
Be  = Gam_f*tAppinv*Hp;

%==========================================================================
%------------------------- the closed-loop system -------------------------
%==========================================================================
% (3a) u_f = Mu*(-Guu u_v - Guy y_v +        r_f - Tcf Hf e_f + Tcf Be e_p)
% (3b) y_f =      Gyu u_v + Gyy y_v + Tuf Mu r_f +  My Hf e_f -  My Be e_p
% where
%   v = varrho = max(rho,p)
%   u_v = past inputs               (nu*v | 1)
%   y_v = past outputs              (ny*v | 1)
%   Mu  = (I + Tcf Tuf)^{-1}
%   My  = (I + Tuf Tcu)^{-1} = I - Tuf Mu Tcf
%   Guu s.t. Guu u_v =         Lcur u_rho + Tcf Lup u_p
%   Guy s.t. Guy y_v =         Lcyr y_rho + Tcf Lyp y_p
%   Gyu s.t. Gyu u_v = -Tuf Mu Lcur u_rho +  My Lup u_p
%   Gyy s.t. Gyy y_v = -Tuf Mu Lcyr y_rho +  My Lyp y_p

% defining matrices
invMu = eye(nu*f) + Tcf*Tuf;  % = I + Tcf Tuf
Mu = inv(invMu);
My = eye(ny*f) - Tuf/invMu*Tcf;
Guu = [zeros(f*nu,(varrho-p)*nu)     Tcf*Lup] + ...
      [zeros(f*nu,(varrho-rho)*nu)   Lcur];
Guy = [zeros(f*nu,(varrho-p)*ny)     Tcf*Lyp] + ...
      [zeros(f*nu,(varrho-rho)*ny)   Lcyr];
Gyu = [zeros(f*ny,(varrho-p)*nu)     My*Lup] + ...
      [zeros(f*ny,(varrho-rho)*nu)  -Tuf/invMu*Lcur];
Gyy = [zeros(f*ny,(varrho-p)*ny)     My*Lyp] + ...
      [zeros(f*ny,(varrho-rho)*ny)  -Tuf/invMu*Lcyr];

xCz_iv_actual = make_ext_ss(nu,ny,f,varrho,invMu,Mu,Tuf/invMu,Guu,Guy,Gyu,Gyy);

%% ===================== estimating system matrices =======================
% steps:
% S1) CL-regression: regress past IO and Rf on [Uf;Yf]. See eqs (3) above
% S2) get Mu & Tuf*Mu from S1
% S3) estimate Tuf: Tuf = (Tuf*Mu)/Mu
% S4) recompute Mu: Mu = (I + Tcf Tuf)^{-1}
% S5) estimate Lup, Lyp
% S6) construct estimates of Gyu, Gyy, Guu, Guy

% make hankel matrices
[~,Uv,Uf] = make_Hankel(u,varrho,f);    
[~,Yv,Yf] = make_Hankel(y,varrho,f);
[~,Rv,Rf] = make_Hankel(wr,varrho,f);
[~,Ev,Ef] = make_Hankel(e,varrho,f);
Ncol = size(Uv,2);

% S1) ------------------------ CL-regression ------------------------------
Wv = [Uv;Yv];  Wp = [Uv(end-nu*p+1:end,:);Yv(end-nu*p+1:end,:)];
WvRf = [Wv;Rf];
MuILuy_est = [Uf;Yf]*WvRf.'/(WvRf*WvRf.');
MuLu_est = MuILuy_est(1:nu*f,:);
Ly_est = MuILuy_est(nu*f+1:end,:);

% S2a) ----------------------- get Mu -------------------------------------
Mu_est1 = MuLu_est(:,end-nu*f+1:end);
Mu_est1 = Mu_est1.*kron(tril(ones(f)),ones(nu)); % impose causal structure
Mu_est1 = blk_toeplitz_mean(Mu_est1,nu,nu);      % avg. over blk-diags
invMu_est1 = inv(Mu_est1);

% S2b) ----------------------- get Tuf*Mu ---------------------------------
TuMu_est1 = Ly_est(:,end-nu*f+1:end);
TuMu_est1 = TuMu_est1.*kron(tril(ones(f)),ones(ny,nu)); % impose causal structure
TuMu_est1 = blk_toeplitz_mean(TuMu_est1,ny,nu);         % avg. over blk-diags

% S3) ------------------------ get Tuf ------------------------------------
Tu_est2 = TuMu_est1/Mu_est1;

% S4) ------------------------ get Mu -------------------------------------
invMu_est2 = eye(f*nu)+Tcf*Tu_est2;
Mu_est2 = inv(invMu_est2);

% S5) ------------------------ get Lup, Lyp -------------------------------
TuMu_est2 = Tu_est2/invMu_est2;
My_est = eye(f*ny) - TuMu_est2*Tcf;

% uf = rf - Lcu_l ul - Lcy_l yl - Tc*yf
% yf = Lu up + Ly yp + Tu uf + H ef + B ep
% rf - Lclu_l ul - Lcy_l yl - (I + Tc Tu) uf = Tc (Lu up + Ly yp + H ef + B ep)
% |__________________ LHS _________________|           \ RHS /

% estimating [Tc;My]*[Lup Lyp] -> Lup & Lyp
UY_v_part = Lcur*Uv(end-nu*rho+1:end,:) + Lcyr*Yv(end-nu*rho+1:end,:);

LHS_u =            Rf -           UY_v_part - invMu_est2*Uf;
LHS_y = -TuMu_est2*Rf + TuMu_est2*UY_v_part + Yf;
TcMyLuy_est = [LHS_u;LHS_y]*Wp.'/(Wp*Wp.');
Luy_est = pinv([Tcf;My_est])*TcMyLuy_est;

Lup_est = Luy_est(:,1:nu*p);
Lyp_est = Luy_est(:,nu*p+(1:ny*f));

% S6) ------------------ get Gyu, Gyy, Guu, Guy ---------------------------
Guu_est = [zeros(f*nu,(varrho-p)*nu)     Tcf*Lup_est] + ...
          [zeros(f*nu,(varrho-rho)*nu)   Lcur];
Guy_est = [zeros(f*nu,(varrho-p)*ny)     Tcf*Lyp_est] + ...
          [zeros(f*nu,(varrho-rho)*ny)   Lcyr];

Gyu_est = [zeros(f*ny,(varrho-p)*nu)     My_est*Lup_est] + ...
          [zeros(f*ny,(varrho-rho)*nu)  -TuMu_est2*Lcur];
Gyy_est = [zeros(f*ny,(varrho-p)*ny)     My_est*Lyp_est] + ...
          [zeros(f*ny,(varrho-rho)*ny)  -TuMu_est2*Lcyr];

%% ===================== construct ext. state-space =======================
xCz_iv = make_ext_ss(nu,ny,f,varrho,invMu_est2,Mu_est2,TuMu_est2,Guu_est,Guy_est,Gyu_est,Gyy_est);

%% plotting

figure(1);
tl = tiledlayout(3,1,'TileSpacing','compact');
nexttile;
imagesc_vik(invMu)
nexttile;
imagesc_vik(invMu_est1)
nexttile;
imagesc_vik(invMu_est2)
set(gcf,'Position',[213   478   560   420]);
title(tl,'$I+\mathcal{T}^\mathrm{c}\mathcal{T}^\mathrm{u}$','Interpreter','latex');

figure(2);
tl2 = tiledlayout(3,1,'TileSpacing','compact');
nexttile;
imagesc_vik(Mu)
nexttile;
imagesc_vik(Mu_est1)
nexttile;
imagesc_vik(Mu_est2)
set(gcf,'Position',[775   480   560   420]);
title(tl2,'$(I+\mathcal{T}^\mathrm{c}\mathcal{T}^\mathrm{u})^{-1}$','Interpreter','latex');

figure(3);
tl3 = tiledlayout(3,1,'TileSpacing','compact');
nexttile;
imagesc_vik(Tuf); ylabel('actual');
nexttile;
imagesc_vik(Tu_est2); ylabel('estimate');
nexttile;
Tu2 = Tuf.*(tril(ones(size(Tuf)),0)+triu(nan(size(Tuf)),1));
dTu = (Tu2-Tu_est2)./Tu2*100;
dTu(isinf(dTu)) = nan;
imagesc_vik(dTu); ylabel('difference %')
set(gcf,'Position',[1336   480   560   420]);
title(tl3,'$\mathcal{T}^\mathrm{u}$','Interpreter','latex');

%% method 1: by Hankel matrices

% [U_iv; Y_iv] = L_cl*[Wv;Rf]
% L_cl = [-Mu*Guu, -Mu*Guy,     Mu;
%             Gyu,     Gyy, Tuf*Mu];

L_cl_est = [-invMu_est2\Guu_est, -invMu_est2\Guy_est, Mu_est2;...
                        Gyu_est,             Gyy_est, TuMu_est2];

L_cl_act = [-invMu\Guu, -invMu\Guy, Mu;...
                Gyu,     Gyy, Tuf/invMu];

% ------------------------ get approximate IV -----------------------------

% compute 'noiseless' uf & yf
UfYf_iv = L_cl_act*WvRf;
Uf_iv = UfYf_iv(1:nu*f,:);
Yf_iv = UfYf_iv(nu*f+1:end,:);


%% method 2: by simulation

end

function SS_ext = make_ext_ss(nu,ny,f,varrho,invMu,Mu,TufMu,Guu,Guy,Gyu,Gyy)
Su = [zeros(varrho*nu,f*nu) eye(varrho*nu)];
Sy = [zeros(varrho*nu,f*nu) eye(varrho*nu)];

Ce_uf = [-invMu\Guu, -invMu\Guy]; De_uf = Mu;
Ce_yf = [Gyu, Gyy];               De_yf = TufMu;
Ce = [Ce_uf;Ce_yf];               De = [De_uf; De_yf];

AeBe_uf = Su*[eye(varrho*nu,varrho*(nu+ny)+nu*f);
              Ce_uf De_uf];
AeBe_yf = Sy*[circshift( eye(varrho*ny,varrho*(nu+ny)+nu*f), varrho*nu, 2 );
              Ce_yf De_yf];
AeBe = [AeBe_uf;AeBe_yf];

Ae = AeBe(:,1:varrho*(nu+ny));
Be = AeBe(:,end-nu*f+1:end);

SS_ext = ss(Ae,Be,Ce,De,[]);
end