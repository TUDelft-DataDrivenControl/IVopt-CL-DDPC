function [opts, u0, y0, xcl0, ...
          Z, Lf, Cz, Tcl, u_cl, y_cl, Cases, u_iv, y_iv, ...
          cost_u1, cost_u2, cost_u, cost_y, cost_tot, FroIDerror] ... 
          = run_sims(opts,sigs,plant,Cz0,Tcl0,yr0,e0,yr1,ur1,e1)
%% This function will run all necessary simulations given input parameters
[A,~,~,D,~] = plant2ABCDK(plant);
[ny,nu] = size(D);
nx = size(A,1);
[p,f,N,Ncl, dRk, Rk, Qk] = deal(opts.p,opts.f,opts.N,opts.Ncl, opts.dRk, opts.Rk, opts.Qk);

% (CL-)SPC weights
dRk= dRk*eye(nu);  dR = kron(speye(f),dRk);
Rk = Rk*eye(nu);   R = kron(speye(f),Rk);
Qk = Qk*eye(ny);   Q = kron(speye(f),Qk);

%% Initial closed-loop data generation
fprintf('Obtaining initial closed-loop data...\n');

% ----------------------- simulation with noise (e0) ----------------------
[uy0,~,xcl0] = lsim(Tcl0,[yr0;e0],[]); uy0 = uy0.'; xcl0 = xcl0.';
u0 = uy0(1:nu,:); y0 = uy0(nu+1:end,:); clear uy0; % get inputs and outputs
xcl0_plant = xcl0(size(Cz0.A,1)+1:size(Cz0.A,1)+nx,end);  % get final state of plant

% ------------------- simulations without future noise --------------------
Uf_iv2 = nan(nu*f,N);
Yf_iv2 = nan(ny*f,N);
tic
for kN = 1:N
    kk = p+kN;
    uy_f_iv_opt = lsim(Tcl0,[yr0(:,kk:kk+f-1);zeros(ny,f)],[],xcl0(:,kk).').';
    Uf_iv2(:,kN) = reshape(uy_f_iv_opt(1:nu,:),nu*f,1);
    Yf_iv2(:,kN) = reshape(uy_f_iv_opt(nu+(1:ny),:),ny*f,1);
end
toc

%% Instrumental Variable Matrices
opts.rho = ss2lag(Cz0); % determine lag of the initial controller

% ============================ create IVs =================================
[Z,Cases,nlcf]= get_Z(u0,y0,yr0,Cz0,Uf_iv2,Yf_iv2,opts); % Z is IV_4_DDPC object containing IVs
nCz = numel(Cases);

% ============= calc. mean & std. dev. of IV trajectories =================
[u_iv,y_iv] = get_uy_iv(Z,nlcf,yr0,opts);

%% Get Subspace Predictive Controllers
% ---------------------- get (CL-)SPC controllers -------------------------
fprintf('Obtaining controllers...\n');

[Lf,Cz,x1_0_SPC,sigs] = get_Lf_Cz(Z,Cases,plant,u0,y0,ur1,yr1,Q,R,dR,opts,sigs);

%% Run closed-loop simulations
fprintf('Running closed-loop simulations...\n');

% ================= create closed-loop systems ============================
Tcl_in  = [sigs.urkf(:)',sigs.yrkf(:)',sigs.ek(:)'];
Tcl_out = [sigs.uk(:)',sigs.yk(:)'];
conOpts = connectOptions("Simplify",false);
for kCz = 1:nCz
    Czn = Cases{kCz};
    Tcl.(Czn) = connect(Cz.(Czn),plant,Tcl_in,Tcl_out,conOpts);
end

% ================= run closed-loop simulations ===========================
% get initial state of closed-loop system
x1_0_plant = plant.A*xcl0_plant+plant.B*[u0(:,end);e0(:,end)]; % initial state of plant
x1_0_cl = [x1_0_SPC;x1_0_plant];

% run simulations
for kCz = 1:nCz
    Czn = Cases{kCz}; % controller name
    uy_cl = lsim(Tcl.(Czn),[ur1(:,f+1:end);yr1(:,f+1:end);e1],[],x1_0_cl).';
    u_cl.(Czn) = uy_cl(1:nu,:);
    y_cl.(Czn) = uy_cl(nu+1:end,:);
end
fprintf('Closed-loop simulations finished!\n');

%% Analysis 
% =================== Calculate average stage costs =======================
fprintf('Calculate average stage costs...\n');

% -------------------------- create functions -----------------------------
% average of (u_k-ur_k).' * Rk * (u_k-ur_k) over Ncl steps
cost_fun_uk = @(u) reshape(ur1(:,1:Ncl)-u,[],1).'*kron(speye(Ncl),Rk)*reshape(ur1(:,1:Ncl)-u,[],1);

% average of du_k.' * dRk * du_k over Ncl steps
cost_fun_duk = @(u) reshape(u-[u0(:,end) u(:,2:end)],[],1).'*kron(speye(Ncl),dRk)*reshape(u-[u0(:,end) u(:,2:end)],[],1);

% average of (y_k-yr_k).' * Rk * (y_k-yr_k) over Ncl steps
cost_fun_yk = @(y) reshape(yr1(:,1:Ncl)-y,[],1).'*kron(speye(Ncl),Qk)*reshape(yr1(:,1:Ncl)-y,[],1);

% -------------------------- perform calculations -------------------------
for kCz = 1:nCz
    Czn = Cases{kCz};
    
    cost_u1.(Czn) = cost_fun_uk( u_cl.(Czn))/Ncl;
    cost_u2.(Czn) = cost_fun_duk(u_cl.(Czn))/Ncl;
    cost_y.(Czn)  = cost_fun_yk( y_cl.(Czn))/Ncl;
    cost_u.(Czn)  = cost_u1.(Czn) + cost_u2.(Czn);
    cost_tot.(Czn)= cost_u.(Czn)  + cost_y.(Czn);  % total cost
end

% ===================== Identification error analysis =====================
fprintf('Calculate identification errors...\n');

for kCz = 1:nCz
    Czn = Cases{kCz};

    IDerror = Lf.(Czn)-Lf.actLf;
    cols = 1:p*nu;
    FroIDerror.(Czn).Up = norm(IDerror(:,cols)); % part affecting Up
    cols = p*nu+(1:p*ny);
    FroIDerror.(Czn).Yp = norm(IDerror(:,cols)); % part affecting Yp
    cols = p*(nu+ny)+(1:f*nu);
    FroIDerror.(Czn).Uf = norm(IDerror(:,cols)); % part affecting Uf
end

end