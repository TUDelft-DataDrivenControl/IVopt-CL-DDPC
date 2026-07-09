function [Lf,Cz,x1_0_SPC,sigs] = get_Lf_Cz(Z,Cases,plant,u0,y0,ur1,yr1,Q,R,dR,opts,sigs)
[nu,ny,p,f,DMCS] = deal(opts.nu,opts.ny,opts.p,opts.f,opts.DMCS);
nCz = numel(Cases);
Yf_r01 = make_Page(y0(:,p+1:end),f,DMCS);

% make functions to get SPC controllers
usol_funs = get_solver(opts,Q,R,dR); % structure w/ up2usol, yp2usol, urf2usol, yrf2usol

% name channels in SPC controllers
sigs.urkf = arrayfun(@(j) sprintf('ur%d_%d',f, j), 1:nu, 'UniformOutput', false);
sigs.yrkf = arrayfun(@(j) sprintf('yr%d_%d',f, j), 1:ny, 'UniformOutput', false);

% for initial state of SPC controllers
up0 = u0(:,end-p+1:end); up0 = up0(:);  % past input data
yp0 = y0(:,end-p+1:end); yp0 = yp0(:);  % past output data
urf = ur1(:,1:f);        urf = urf(:);  % future input references
yrf = yr1(:,1:f);        yrf = yrf(:);  % future output references

% create controllers
for iCz = 1:nCz
    Czn = Cases{iCz}; % name of controller

    % ============================ get Lf (estimate) ==========================
    % -------------------------- (1-6) SPCs based on an IV --------------------
    if startsWith(Czn,'iv')
        Lf.(Czn) = Yf_r01*Z.(Czn).'*pinv(Z.iv1*Z.(Czn).');

    % -------------------------- (7) CL-SPC -----------------------------------
    elseif strcmp(Czn, 'CLSPC')
        Lf.(Czn) = get_Lf_CL_SPC(u0,y0,p,f,nu,ny,DMCS);

    % -------------------------- (8) SPC w/ actual Lf -------------------------
    elseif strcmp(Czn, 'actLf')
        [Up2Yf_inno, Yp2Yf_inno, Uf2Yf_inno] = get_actual_matrices(plant,p,f);
        Lf.(Czn) = [Up2Yf_inno Yp2Yf_inno Uf2Yf_inno];

    % -------------------------- (9) Transient Predictor ----------------------
    elseif strcmp(Czn, 'TrPred')
        Lf.(Czn) = get_Lf_TransPred(u0,y0,p,f,nu,ny,DMCS);
    end

    % =========================== get controller ==============================
    if iCz > 1
        Cz.(Czn) = Lf_2_SPC(Lf.(Czn),usol_funs,opts,sigs);
    else
        % also create initial state of the controller
        [Cz.(Czn),x1_0_SPC] = Lf_2_SPC(Lf.(Czn),usol_funs,opts,sigs,up=up0,yp=yp0,urf=urf,yrf=yrf);
    end
end
end

function Lf = get_Lf_CL_SPC(u1,y1,p,f,nu,ny,DMCS)
    Upf = make_Page(u1,p+f,DMCS); Up = Upf(1:nu*p,:); Uf = Upf(nu*p+1:nu*(p+1),:);
    Ypf = make_Page(y1,p+f,DMCS); Yp = Ypf(1:ny*p,:); Yf = Ypf(ny*p+1:ny*(p+1),:);
    L1 = Yf*pinv([Up;Yp;Uf(1:nu,:)]);
    
    C_tKpu_hat = L1(:,1:p*nu);
    C_tKpy_hat = L1(:,p*nu+1:p*(nu+ny));
    D_hat = L1(:,end-nu+1:end);
    
    % construct multiple-step-ahead predictor
    tLest_u = zeros(ny*f,nu*(p+f));
    tLest_u(1:ny,1:nu*(p+1)) = [C_tKpu_hat D_hat];
    tLest_y = zeros(ny*f,ny*(p+f));
    tLest_y(1:ny,1:ny*p) = C_tKpy_hat;
    for kr = 2:f
        tLest_u((kr-1)*ny+1:kr*ny,:) = circshift(tLest_u((kr-2)*ny+1:(kr-1)*ny,:),nu,2);
        tLest_y((kr-1)*ny+1:kr*ny,:) = circshift(tLest_y((kr-2)*ny+1:(kr-1)*ny,:),ny,2);
    end
    tHf = eye(ny*f)-tLest_y(:,end-ny*f+1:end);
    Lf = tHf\[tLest_u(:,1:p*nu) tLest_y(:,1:p*ny) tLest_u(:,end-nu*f+1:end)];
end

function Lf = get_Lf_TransPred(u1,y1,p,f,nu,ny,DMCS)
    Zpf = make_Page([u1;y1],p+f,DMCS); Zp = Zpf(1:(nu+ny)*p,:); Zf = Zpf((nu+ny)*p+1:end,:);

    [~,R] = qr([Zp;Zf].','econ'); R = R.';
    tLest = zeros(ny*f,(p+f)*(nu+ny));
    for kf = 1:f
        % get R11 & R21
        Rzp_idx = 1:(p+kf-1)*(nu+ny)+nu;
        Ryk_idx = Rzp_idx(end)+(1:ny);
        R11 = R(Rzp_idx,Rzp_idx);
        R21 = R(Ryk_idx,Rzp_idx);

        % compute relevant part of tLest matrix
        rows = (kf-1)*ny+1:kf*ny;
        cols = 1:p*(nu+ny)+kf*nu+(kf-1)*ny;
        tLest(rows,cols) =  R21*pinv(R11);
    end
    ucols = mod(1:size(tLest,2),nu+ny); ucols(ucols == 0) = nu+ny;
    msk_u = ucols < nu+1;
    msk_y = ~msk_u;
    tLest_u = tLest(:,msk_u);
    tLest_y = tLest(:,msk_y);
    tHf = eye(ny*f)-tLest_y(:,end-ny*f+1:end);
    Lf = tHf\[tLest_u(:,1:p*nu) tLest_y(:,1:p*ny) tLest_u(:,end-nu*f+1:end)];
end