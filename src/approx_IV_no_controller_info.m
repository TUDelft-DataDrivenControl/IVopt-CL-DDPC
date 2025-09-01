function [Uf_iv,Yf_iv] = approx_IV_no_controller_info(u,y,w,rho,opts)
% This function approximates the IV matrices Uf and Yf that would have
% been obtained without future noise Ef, using no knowledge of the
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

[p,f,nu,ny] = deal(opts.p,opts.f,opts.nu,opts.ny);
varrho = max(rho,p);

[~,Nbar] = size(u);
validateattributes(u, {'double'},{'size',[nu Nbar]});
validateattributes(y, {'double'},{'size',[ny Nbar]});
validateattributes(w, {'double'},{'size',[ny Nbar]});

% create data matrix D and Uf, Yf
[~,Uv,Uf] = make_Hankel(u,varrho,f);
[~,Yv,Yf] = make_Hankel(y,varrho,f);
[~,Wr,Wf] = make_Hankel(w,varrho,f);
Wr = Wr(end-rho+1:end,:); % select last rho block rows
D = [Uv;Yv;Wr;Wf];

%% ------------------------- get closed-loop system ------------------------
L_cl = [Uf;Yf]*pinv(D); % initial estimate

% split up into parts
L_clu = L_cl(1:nu*f,:);     % D -> uf
L_cly = L_cl(nu*f+1:end,:); % D -> yf

% modify parts representing influence of Wf (impose causality & average)
L_clur = L_clu(:,end-f*ny+1:end);
L_clyr = L_cly(:,end-f*ny+1:end);

% part for Wf should be causal -> impose (block-)tril structure
L_clur = L_clur.*kron(tril(ones(f)),ones(nu,ny)); % blocks: nu x ny
L_clyr = L_clyr.*kron(tril(ones(f)),ones(ny,ny)); % blocks: ny x ny

% take mean of corresponding elements in blocks along diagonals
L_clur = blk_toeplitz_mean(L_clur,nu,ny);
L_clyr = blk_toeplitz_mean(L_clyr,ny,ny);

% adjust L_clu & L_cly accordingly
L_clu(:,end-f*ny+1:end) = L_clur;
L_cly(:,end-f*ny+1:end) = L_clyr;

% concatenate averaged matrices
L_cl = [L_clu; L_cly];

%% ------------------------ get approximate IV -----------------------------
UfYf_iv = L_cl*D;
Uf_iv = UfYf_iv(1:nu*f,:);
Yf_iv = UfYf_iv(nu*f+1:end,:);
% [Rpf2,~,~] = make_Hankel([zeros(ny,s) ref],s,pf); % initialize by zero
% 
% % compute 'noiseless' u & y
% UY_iv = L_cl*Rpf2;
% U_iv = UY_iv(1:nu*pf,:);
% Y_iv = UY_iv(nu*pf+1:end,:);
% 
% % average over found u & y
% U_iv = flipud(blk_toeplitz_mean(flipud(U_iv),nu,1));
% Y_iv = flipud(blk_toeplitz_mean(flipud(Y_iv),ny,1));
% 
% % select u & y from block-Hankel matrices as IVs
% u_iv = U_iv(:,1); u_iv = [reshape(u_iv,nu,[]),U_iv(end-nu+1:end,2:end)];
% y_iv = Y_iv(:,1); y_iv = [reshape(y_iv,ny,[]),Y_iv(end-ny+1:end,2:end)];
end