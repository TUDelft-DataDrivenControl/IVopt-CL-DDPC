function [u_iv,y_iv,xCz_iv] = approx_IV_controller_info(u,y,ref,p,f,Ziv_init,Cz,nit_max)
[nu,Nbar] = size(u);
ny = size(y,1);
validateattributes(y,  {'double'},{'size',[ny Nbar]});
validateattributes(ref,{'double'},{'size',[ny Nbar]});

[~,Up,Uf] = make_Hankel(u,p,f);
[~,Yp,Yf] = make_Hankel(y,p,f);

ref = [ref zeros(ny,ceil(Nbar/ny/f)*ny*f-Nbar)]; % padding ref with zeros for simulation
Ref = reshape(ref,ny*f,[]);
Wp  = [Up;Yp;Uf];
Ziv = Ziv_init;

%% construct constant matrices
% controller-dependent matrices
[Ac,Bc,Cc,Dc] = ssdata(Cz);
Tf_c = make_blk_tril_toeplitz(Ac,Bc,Cc,Dc,f);
Gf_c = make_ext_obsv(Ac,Cc,f);
Kf_c = make_ext_ctrb(Ac,Bc,f,rev=true);
Acf  = Ac^f;

% get predictor form with eigs(Ac-Kc*Cc) \approx 0
nxc = size(Ac,1);
Kc = acker(Ac.',Cc.',zeros(nxc,1)).';
tAc = Ac-Kc*Cc;
tBc = Bc-Kc*Dc;
tKp_ce = make_ext_ctrb(tAc,tBc,p,rev=true);
tKp_cu = make_ext_ctrb(tAc,Kc,p,rev=true);

% selection matrices
Su = [zeros(p*nu,f*nu) eye(p*nu)];
Su1 = Su(:,1:p*nu);     % nu*p x nu*p
Su2 = Su(:,p*nu+1:end); % nu*p x nu*f
Sy = [zeros(p*ny,f*ny) eye(p*ny)];
Sy1 = Sy(:,1:p*ny);     % ny*p x ny*p
Sy2 = Sy(:,p*ny+1:end); % ny*p x ny*f

% actual system:
% [plant,nx,nu,ny,A,B,C,D,K,~] = model_Landau1995(Re=4.81e-1);
% tA = A-K*C;
% tB = B-K*D;
% Obsv_plant = make_ext_obsv(A,C,f);
% alpha_u = make_ext_ctrb(tA,tB,p,rev=true);-tA^p*pinv(make_ext_obsv(A,C,p))*make_blk_tril_toeplitz(A,B,C,D,p);
% alpha_u = Obsv_plant*alpha_u;
% alpha_y = make_ext_ctrb(tA,K,p,rev=true);+tA^p*pinv(make_ext_obsv(A,C,p));
% alpha_y = Obsv_plant*alpha_y;
% Tf_u = make_blk_tril_toeplitz(A,B,C,D,f);

I1 = eye(nu*f);
I2 = eye(ny*f);
%% Iterations
for kit = 1:nit_max

% ------------------------------ identification ---------------------------
% Lest = Yf*Ziv.'*pinv([Up;Yp;Uf]*Ziv.');
Ziv = Wp*Ziv.'/(Ziv*Ziv.')*Ziv;
Lest = Yf*Ziv.'/(Ziv*Ziv.');
% % Lest = Yf*pinv([Up;Yp;Uf]);
alpha_u = Lest(:,1:p*nu);
alpha_y = Lest(:,p*nu+1:p*(nu+ny));
Tf_u = Lest(:,end-nu*f+1:end);
Tf_u = blk_toeplitz_mean(Tf_u,ny,nu);
% Tf_u = Tf_u.*kron(tril(ones(f),-3),ones(ny,nu));

% ---------------------- create closed-loop system ------------------------
% x_j = [ u_{k-p,k-1};
%         y_{k-p,k-1};
%         x_k^c];
lam1 = I2+Tf_u*Tf_c;
E11 = [eye(nu*f) Tf_c;...
       -Tf_u eye(ny*f)];
invE11 = [I1-Tf_c/lam1*Tf_u -Tf_c/lam1;
                  lam1\Tf_u   inv(lam1)];
A11 = [Gf_c*tKp_cu -Gf_c*tKp_ce;...
       alpha_u alpha_y];

if any(abs(eig(A11)) > 1)
    [Veigs,Deigs] = eig(A11);
    eigs = diag(Deigs);
    msk = abs(eigs) > 1;    
    eigs(msk) = conj(1./eigs(msk));

    Deigs_new = diag(eigs);
    A11_new = real(Veigs*Deigs_new/Veigs);
    A11(end-ny*f+1:end,:) = A11_new(end-ny*f+1:end,:);
end
Gf_c_tKp_ce = -A11(1:ny*f,nu*p+1:(nu+ny)*p);
A = invE11*[A11 [Gf_c_tKp_ce; zeros(nu*f,ny*p)]];
C = A;
A = [A;zeros(ny*f,p*(nu+2*ny))];


% invEA11 = invE11*A11;
% A = blkdiag(invEA11,zeros(ny*f));


% A(1:ny*f,end-p*ny+1:end) = -A(1:ny*f,nu*p+1:(nu+ny)*p);%Gf_c*tKp_ce;
% C = A(1:(nu+ny)*f,:);

% define useful matrices
lam1  = eye(f*nu)+Tf_c*Tf_u;
lam5  = lam1\Tf_c;%Tf_c/(eye(f*nu)+Tf_u*Tf_c);%
lam2  = -lam5*alpha_u;
lam3  = -lam5*alpha_y;
lam4  = -lam1\Gf_c;
lam6  = alpha_u + Tf_u*lam2;
lam7  = alpha_y + Tf_u*lam3;
lam8  = Tf_u*lam4;
lam9  = Tf_u*lam5;

% create state-space
% A = [ Su1+Su2*lam2       Su2*lam3        Su2*lam4;...
%           Sy2*lam6   Sy1+Sy2*lam7        Sy2*lam8;...
%         -Kf_c*lam6     -Kf_c*lam7   Acf-Kf_c*lam8];
% B = [       Su2*lam5;...
%             Sy2*lam9;...
%       Kf_c-Kf_c*lam9];
% C = [         lam2           lam3            lam4;...
%               lam6           lam7            lam8];
% D = [           lam5;...
%                 lam9];

% Ta = blkdiag(eye(p*nu),eye(p*ny),tKp_ce);
% Ta(end-nxc+1:end,:) = [tKp_cu -tKp_ce tKp_ce];
% 
% A0 = [ Su1+Su2*lam2       Su2*lam3        Su2*lam4;...
%           Sy2*lam6   Sy1+Sy2*lam7        Sy2*lam8];
% A = [A0*Ta; zeros(ny*p,p*(nu+ny)) Sy1];

B = [       Su2*lam5;...
            Sy2*lam9;...
            Sy2];
% C = [         lam2           lam3            lam4;...
%               lam6           lam7            lam8];
% C = C*Ta;
D = [           lam5;...
                lam9];

% A(1:p*(nu+ny),1:p*nu-nxc+1) = 0;
% C(:,1:p*nu-nxc+1) = 0;
% A(1:p*(nu+ny),p*nu+1:p*(nu+ny)-nxc+1) = 0;
% C(:,p*nu+1:p*(nu+ny)-nxc+1) = 0;
% A(1:p*(nu+ny),p*(nu+ny)+1:end-nxc+1) = 0;
% C(:,p*(nu+ny)+1:end-nxc+1) = 0;

sys_cl = ss(A,B,C,D,[]);

% initial state
x0_cl = zeros(size(A,1),1);%[zeros(nu*p,1);zeros(ny*p,1);zeros(ny*p,1)];

% ------------- ensure stability of closed-loop system --------------------
eigs = pole(sys_cl);
itk = 0;
while any(abs(eigs) > 1)
    % analysis A matrix of closed-loop system
    % plotABCD2(A,B,C,D,p,f,nu,ny,nxc);

    itk = itk + 1;

    A11 = A(1:p*(nu+ny),1:p*(nu+ny));
    [Veigs,Deigs] = eig(A11);
    eigs = diag(Deigs);
    msk = abs(eigs) > 1;    
    eigs(msk) = conj(1./eigs(msk));

    Deigs_new = diag(eigs);
    A11_new = Veigs*Deigs_new/Veigs;
    A(1:p*(nu+ny),1:p*(nu+ny)) = real(A11_new);
end

% adjust C matrix if A has been adjusted
if itk > 0
    % determine how many last block-rows are affected
    nrows = min(f,p);
    C(nu*f-nrows*nu+1:nu*f,:) = A(nu*p-nrows*nu+1:nu*p,:);
    C(end-nrows*ny+1:end,:)   = A((nu+ny)*p-nrows*ny+1:(nu+ny)*p,:);
    sys_cl = ss(A,B,C,D,-1);

    % plotABCD2(A,B,C,D,p,f,nu,ny,nxc)
end

% ------------------- get approximate optimal IV --------------------------
% uin = reshape(u(:,1:floor(Nbar/f)*f),nu*f,[]);
% yin = reshape(y(:,1:floor(Nbar/f)*f),ny*f,[]);
% cl_data = iddata([uin;yin].',Ref(:,1:end-1).',1);
% x0_cl = findstates(idss(sys_cl),cl_data);
% simulate closed-loop system
[uy_iv,~,x_iv] = lsim(sys_cl,Ref.',[],x0_cl); uy_iv = uy_iv.'; x_iv = x_iv.';
xCz_iv = x_iv(end-size(Cz.A,1)+1:end,:);
% x_iv_new = sys_cl.A*x_iv(:,end) + sys_cl.B*Ref(:,end);
% uy_iv = [uy_iv sys_cl.C*x_iv_new+sys_cl.D*Ref(:,end)];

% get approximate optimal IV
u_iv = reshape(uy_iv(1:f*nu,:),nu,[]);
y_iv = reshape(uy_iv(f*nu+1:end,:),ny,[]);

u_iv = u_iv(:,1:Nbar); % limit to last Nbar values (reference was padded)
y_iv = y_iv(:,1:Nbar);

if kit < nit_max
    % construct IV for further iterations
    [~,Up_iv,Uf_iv] = make_Hankel(u_iv,p,f);
    [~,Yp_iv,~] = make_Hankel(y_iv,p,f);
    Ziv = [Up_iv;Yp_iv;Uf_iv];
end

end
end

function plotABCD2(A,B,C,D,p,f,nu,ny,nxc)
fig5 = figure(5); tl5 = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
ax5_1 = nexttile(tl5);
spy([A,B;C,D]); grid on;
xyline(p*nu+0.5);             % up
xyline(p*(nu+ny)+0.5);        % yp
xyline2(p*(nu+ny)+nxc+0.5,2); % rp
xyline(p*(nu+ny)+nxc+f*nu);
ax5_2 = nexttile(tl5);
imagesc_vik([A,B;C,D]);
xyline(p*nu+0.5);             % up
xyline(p*(nu+ny)+0.5);        % yp
xyline2(p*(nu+ny)+nxc+0.5,2); % rp
xyline(p*(nu+ny)+nxc+f*nu);
fig5.Position = [220 151 904 419];

    function xyline(x)
        xline(x); yline(x);
    end
    function xyline2(x,width)
        xline(x,'LineWidth',width); yline(x,'LineWidth',width);
    end
end