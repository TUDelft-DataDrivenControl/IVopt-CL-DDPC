%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           DDPC using an Optimal-IV
%           Authors: R. Dinkla, T. Oomen, J.W. van Wingerden
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solves a constrained least-squares regression problem
clear;
close all;
rng default;

%% problem parameters
f = 10;
fid = f;
p = 20;
ny = 2;
nu = 2;
N = 1000;
Nbar = p + fid + N -1;

%% make optimization variables and parameters

% parameters
y_ = casadi.SX.sym('y',ny,Nbar);
u_ = casadi.SX.sym('u',nu,Nbar);
p_ = [u_(:);y_(:)];

% optimization variables
e_ = casadi.SX.sym('e',ny,Nbar);
tCAiB_D_ = casadi.SX.sym('tCAiB_D',ny,(p+f)*nu); % predictor Markov params - u
tCAiK_ = casadi.SX.sym('tCAiK',ny,(p+f-1)*ny);   % predictor Markov params - y
eF_ = e_(:,end-(fid+N-1)+1:end);
x_ = [tCAiB_D_(:);tCAiK_(:);eF_(:)]; % ny*[fid+N-1-ny + (p+f)*(nu+ny)] unknowns

% number of optimization variables in each
neF = (fid+N-1)*ny;
ntCAiB_D = (p+f)*nu*ny;
ntCAiK = ny*(p+f-1)*ny;

% derived variables - data
Ypf_ = make_CasADi_Hankel(y_,p+fid,N,'y');
Yp_ = Ypf_(1:p*ny,:);
Yf_ = Ypf_(p*ny+1:end,:);
Upf_ = make_CasADi_Hankel(u_,p+fid,N,'u');
Up_ = Upf_(1:p*nu,:);
Uf_ = Upf_(p*nu+1:end,:);
Epf_ = make_CasADi_Hankel(e_,p+fid,N,'e');
Ep_ = Epf_(1:p*ny,:);
Ef_ = Epf_(p*ny+1:end,:);

% derived variables - input predictor Markov parameters
tGK_Hu_ = casadi.SX.zeros(f*ny,(p+f)*nu);
for k = f:-1:1
    tGK_Hu_((k-1)*ny+1:k*ny,1:(p+k)*nu) = tCAiB_D_(:,(f-k)*nu+1:end);
end
tGKu_ = tGK_Hu_(:,1:p*nu);
tHu_ = tGK_Hu_(:,end-nu*f+1:end);

% get_GKu = casadi.Function('get_GKu',{tCAiB_D_},{tGKu_});
% CAiB_D = rand(size(CAiB_D_));
% GKu = full(get_GKu(CAiB_D));
% figure();
% imagesc(GKu)
% colorbar

% derived variables - output predictor Markov parameters
tGK_Hy_ = casadi.SX.zeros(f*ny,(p+f)*ny);
tCAiK_0_ = [tCAiK_ casadi.SX.zeros(ny,ny)];
for k = f:-1:1
    tGK_Hy_((k-1)*ny+1:k*ny,1:(p+k)*ny) = tCAiK_0_(:,(f-k)*ny+1:(p+f)*ny);
end
tGKy_ = tGK_Hy_(:,1:p*ny);
tHy_ = tGK_Hy_(:,end-ny*f+1:end);

% get_GKy = casadi.Function('get_GKy',{tCAiK_},{tGKy_});
% CAiK = rand(size(CAiK_));
% GKy = full(get_GKy(CAiK));
% figure();
% imagesc(GKy)
% colorbar

%% Run simulation
% make nominal system
nx = 2;
sys = drss(nx,ny,nu);
[A,B,C,D] = ssdata(sys);

% set noise properties
Re = 1e1*tril(rand(ny)*2); Re = Re*Re.';% 1e2*eye(ny);
K = place(A.',C.',linspace(0.9,0.95,nx)).';

% create system with noise
sys2 = ss(A,[B K],C,[D eye(ny)],-1);

At = A-K*C;
Bt = B-K*D;

tCAiB = cell(1,p+f-1);
tAiB = Bt;
for k = p+f-1:-1:1
    tCAiB{k} = C*tAiB;
    tAiB = At*tAiB;
end
tCAiB = cell2mat(tCAiB);
tCAiB_D = [tCAiB D];

tCAiK = cell(1,p+f-1);
tAiK = K;
for k = p+f-1:-1:1
    tCAiK{k} = C*tAiK;
    tAiK = At*tAiK;
end
tCAiK = cell2mat(tCAiK);

% ------------------------ run simulation ---------------------------
Ru = eye(nu)*1e2; % set OL-input variance
u = mvnrnd(zeros(nu,1),Ru,Nbar).';
e = mvnrnd(zeros(ny,1),Re,Nbar).'; eF = e(:,end-(fid+N-1)+1:end);
y = lsim(sys2,[u;e],[]).';

% make Hankel matrices
get_Ypf = casadi.Function('get_Ypf',{y_},{Ypf_});
Ypf = full(get_Ypf(y));
Yf = Ypf(p*ny+1:end,:);
get_Upf = casadi.Function('get_Upf',{u_},{Upf_});
Upf = full(get_Upf(u));

%% Equation:
% alpha: augmented matrix with unknowns
alpha = [tGK_Hu_ tGK_Hy_ Ef_];
alpha_x = jacobian(alpha(:),x_);
alpha_x = sparse(casadi.DM(alpha_x));

% beta: matrix with all regression parameters
% beta = [Upf_;Ypf_;casadi.DM.eye(N)];
% Ax = kron(beta.',casadi.DM.eye(ny*f))*alpha_x;
% get_Ax = casadi.Function('get_Ax',{p_},{Ax});
% p_test = rand(size(p_));
% Yf = full(get_Yf(p_test));
% Ax_test = full(get_Ax(p_test));
beta = [Upf;Ypf;speye(N)];
Ax = kron(beta.',speye(ny*f))*alpha_x;

%{
% cost
cost = (Yf(:)-Ax*x_).'*(Yf(:)-Ax*x_);
prob    = struct('f', cost, 'x', x_);%,'p',p_);
solver_opts = struct('solver', 'osqp','options', struct('print_time',0));
QP1 = casadi.qpsol( 'solver', solver_opts.solver,prob,solver_opts.options);
res = QP1();
%}

[U,S,V] = svd(full(Ax));
S = sparse(S);
rS = sum(full(diag(S)>max(size(S))*eps(diag(S).'*diag(S)))); % determine the rank
U1 = U(:,1:rS);
U2 = U(:,rS+1:end);
V1 = V(:,1:rS);
V2 = V(:,rS+1:end);
S1 = sparse(S(1:rS,1:rS));

Qw = kron(speye(neF/ny),eye(ny));
V1e = V1(end-neF+1:end,:);
V2e = V2(end-neF+1:end,:);
H = V2e.'*Qw*V2e/2; H = (H+H.')/2;
ct = 2*Yf(:).'*U1/S1*V1e.'*Qw*V2e;

wopt = -H\ct.';%-4*S1\U1.'*Yf(:);
xopt = V1/S1*U1.'*Yf(:) + V2*wopt;

tCAiB_D_est = xopt(1:ntCAiB_D,1);               tCAiB_D_est = reshape(tCAiB_D_est,ny,[]);
tCAiK_est = xopt(ntCAiB_D+1:ntCAiB_D+ntCAiK,1); tCAiK_est   = reshape(tCAiK_est,ny,[]);
eF_est = xopt(end-neF+1:end,1);                 eF_est      = reshape(eF_est,ny,[]);
%%
figure();
t = tiledlayout(3,2,'TileSpacing','compact','Padding','compact');
nexttile;
imagesc(tCAiB_D)
ylabel('$C\tilde{A}^i\tilde{B}$ $D$','interpreter','latex');
nexttile;
imagesc(tCAiB_D_est)

nexttile;
imagesc(tCAiK)
ylabel('$C\tilde{A}^iK$','interpreter','latex');
colorbar
nexttile;
imagesc(tCAiK_est);

nexttile;
imagesc(eF);
ylabel('noise','interpreter','latex');
xlabel('real')
nexttile;
imagesc(eF_est);
xlabel('estimated');


figure(2);
subplot(2,1,1);
plot(eF(1,:)); hold on;
plot(eF_est(1,:));
subplot(2,1,2);
plot(eF(2,:)); hold on;
plot(eF_est(2,:));
legend('real','est')
%% Helper functions
function Hankel = make_CasADi_Hankel(CasADi_var,dim1,dim2,u_or_y)
% construct RHS for equality governing dynamics
var_dim1 = size(CasADi_var,1);

% -> construct RHS block-row by block-row
X = casadi.MX.sym('X',var_dim1,dim1+dim2-1);
ys = cell(dim1,1);
for i=1:dim1
    ys{i,1} = X(:,i:i+dim2-1);
end
Y = vertcat(ys{:});
F = casadi.Function(['F',u_or_y],{X},{Y}); % to ensure distinct function names
Hankel = F(CasADi_var);
end