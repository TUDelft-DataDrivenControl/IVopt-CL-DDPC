function [u_iv,y_iv] = approx_IV_no_controller_info(u,y,ref,s,p,f)
pf = p+f; % total size of 'future' window needed
[nu,Nbar] = size(u);
ny = size(y,1);
validateattributes(y,  {'double'},{'size',[ny Nbar]});
validateattributes(ref,{'double'},{'size',[ny Nbar]});

[~,~,Upf] = make_Hankel(u,s,pf);
[~,~,Ypf] = make_Hankel(y,s,pf);
[Rs,~,~] = make_Hankel(ref,s,pf);

% ------------------------- get closed-loop system ------------------------
L_cl = [Upf;Ypf]*pinv(Rs); % initial estimate

% split up into parts
L_clu = L_cl(1:nu*pf,:);     % rf -> uf
L_cly = L_cl(nu*pf+1:end,:); % rf -> yf

% should be causal -> impose (block-)tril structure
L_clu = L_clu.*kron(tril(ones(pf,s+pf),s),ones(nu,ny)); % blocks: nu x ny
L_cly = L_cly.*kron(tril(ones(pf,s+pf),s),ones(ny,ny)); % blocks: ny x ny

% take mean of corresponding elements in blocks along diagonals
L_clu = blk_toeplitz_mean(L_clu,nu,ny);
L_cly = blk_toeplitz_mean(L_cly,ny,ny);

% concatenate averaged matrices
L_cl = [L_clu; L_cly];
figure(1);
imagesc_vik(L_cl)

% ------------------------ get approximate IV -----------------------------
[Rpf2,~,~] = make_Hankel([zeros(ny,s) ref],s,pf); % initialize by zero

% compute 'noiseless' u & y
UY_iv = L_cl*Rpf2;
U_iv = UY_iv(1:nu*pf,:);
Y_iv = UY_iv(nu*pf+1:end,:);

% average over found u & y
U_iv = flipud(blk_toeplitz_mean(flipud(U_iv),nu,1));
Y_iv = flipud(blk_toeplitz_mean(flipud(Y_iv),ny,1));

% select u & y from block-Hankel matrices as IVs
u_iv = U_iv(:,1); u_iv = [reshape(u_iv,nu,[]),U_iv(end-nu+1:end,2:end)];
y_iv = Y_iv(:,1); y_iv = [reshape(y_iv,ny,[]),Y_iv(end-ny+1:end,2:end)];
end