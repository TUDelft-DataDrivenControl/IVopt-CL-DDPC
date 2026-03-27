function [Uf_iv,Yf_iv] = approx_IV_controller_info(u,y,w,opts,Cz)
% This function approximates the IV matrices Uf and Yf that would have
% been obtained without future noise Ef, using knowledge of the
% functional form of the 1 DOF feedback controller C_{fb}:
% 
% feedback law:
%       u_k = C_{fb}(q) (w_k - y_k)
%
% Resulting composition of data matrix D:
%       D = [U_{\varrho}; Y_{\varrho}; W_{\rho}; W_f]

[rho,p,f] = deal(opts.rho,opts.p,opts.f);
[nu,Nbar] = size(u);
ny = size(y,1);
validateattributes(y, {'double'},{'size',[ny Nbar]});
validateattributes(w, {'double'},{'size',[ny Nbar]});

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
%==========================================================================
%--------------------------- the controller -------------------------------
%==========================================================================
% (1)      u_f = Lce er_rho + Lcu u_rho + Tcf er_f
% where:
%   rho    = lag of the controller
%   u_f    = future inputs                  (nu*f   | 1)
%   u_rho  = past inputs                    (nu*rho | 1)
%   er_rho = w_rho - y_rho (tracking error)
%       w_rho  = past exogenous inputs      (ny*rho | 1)
%       y_rho  = past outputs               (ny*rho | 1)
%   er_f   = w_f - y_f     (tracking error)
%       w_f    = future exogenous inputs    (ny*f   | 1)
%       y_f    = future outputs             (ny*f   | 1)

[Ac,Bc,Cc,Dc] = ssdata(Cz);
[Lce, Lcu, Tcf] = get_lifted_mats(Ac,Bc,Cc,Dc,[],rho,f);

%{
% actually unknown parts:
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

[Lup, Lyp, Tuf, Hf, Be] = get_lifted_mats(A,B,C,D,K,p,f);

%==========================================================================
%------------------------- the closed-loop system -------------------------
%==========================================================================
% (3a) u_f = Guu u_v + Guy y_v + Guw w_vf - Mu Tcf (Hf e_f - Be e_p)
% (3b) y_f = Gyu u_v + Gyy y_v + Gyu w_vf +     My (Hf e_f - Be e_p)
% where
%   v = varrho = max(rho,p)
%   u_v = past inputs               (nu*v | 1)
%   y_v = past outputs              (ny*v | 1)
%   Mu  = (I + Tcf Tuf)^{-1}
%   My  = (I + Tuf Tcu)^{-1} = I - Tuf Mu Tcf
%   Guu s.t. Guu u_v =    Mu ( Lcu u_rho - Tcf Lup u_p )
%   Guy s.t. Guy y_v =    Mu ( Lce y_rho - Tcf Lyp y_p )
%   Gyu s.t. Gyu u_v = My Lup u_p + Tuf Mu Lcu u_rho
%   Gyy s.t. Gyy y_v = My Lyp y_p - Tuf Mu Lce y_rho

% defining matrices
invMu = eye(nu*f) + Tcf*Tuf;  % = I + Tcf Tuf
Mu = inv(invMu);
TufMu = Tuf/invMu;
My = eye(ny*f) - TufMu*Tcf;

Guu = invMu\add_over( Lcu,-Tcf*Lup);
Guy = invMu\add_over(-Lce,-Tcf*Lyp);
Guw = invMu\[Lce Tcf];
Gu = [Guu Guy Guw];

Gyu = add_over(My*Lup,  TufMu*Lcu);
Gyy = add_over(My*Lyp, -TufMu*Lce);
Gyw = Tuf*Guw;
Gy = [Gyu Gyy Gyw];
%}

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
[~,Wr,Wf] = make_Hankel(w,varrho,f); % Wr -> rho instead of varrho deep

% select specific rho and p long matrices
Wr = Wr(end-rho*ny+1:end,:);
Up = Uv(end-p*nu+1:end,:); Ur = Uv(end-rho*nu+1:end,:);
Yp = Yv(end-p*ny+1:end,:); Yr = Yv(end-rho*ny+1:end,:);

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

%% Helper functions
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
        varargout{1}  = make_blk_tril_toeplitz(A,K,C,eye(nyC),f); % = Hf
    end
    if nargout >= 5
        Hp  = make_blk_tril_toeplitz(A,K,C,eye(nyC),p);
        varargout{2}  = Gamf*tAp_pinv*Hp; % = Be
    end
end

function mat3 = add_over(mat1,mat2)
% adds matrices with necessary prepending of zeros
[s11,s12] = size(mat1);
validateattributes(mat2,  {'double'},{'size',[s11 NaN]});

s22 = size(mat2,2);
if s22 >= s12
    mat3 = mat2;
    mat3(:,end-s12+1:end) = mat3(:,end-s12+1:end) + mat1;
else
    mat3 = mat1;
    mat3(:,end-s22+1:end) = mat3(:,end-s22+1:end) + mat2;
end
end

function Mat = blk_tril_avg(Mat,dim1,dim2)
    Mat = make_blk_tril(Mat, [dim1,dim2]);
    Mat = blk_toeplitz_mean(Mat,dim1,dim2);
end

