function [get_usol,varargout] = get_solver(opts,Q,R,dR)
arguments (Input)
    opts struct
    Q        double {mustBeMatrix}
    R        double {mustBeMatrix}
    dR       double {mustBeMatrix}
end
[nu,ny,p,f] = deal(opts.nu,opts.ny,opts.p,opts.f);
% making controller matrices
uf_ = casadi.SX.sym('uf',nu,f);
urf_ = casadi.SX.sym('urf',nu,f);
yrf_ = casadi.SX.sym('yrf',ny,f);
up_ = casadi.SX.sym('up',nu,p);
yp_ = casadi.SX.sym('yp',ny,p);

Gf_tKp_u_ = casadi.SX.sym('Gf_tKpu',ny*f,nu*p);
Tf_u_     = casadi.SX.sym('Tfu',ny*f,nu*f);
Gf_tKp_y_ = casadi.SX.sym('Gf_tKpu',ny*f,ny*p);
Lest_ = [Gf_tKp_u_ Gf_tKp_y_ Tf_u_];
yf_ = Gf_tKp_u_*up_(:) + Tf_u_*uf_(:) + Gf_tKp_y_*yp_(:);

% make casadi Functions
% get_Gf_tKp_u_ = casadi.Function('get_Gf_tKp_u',{Lest_},{Gf_tKp_u_});
% get_Gf_tKp_y_ = casadi.Function('get_Gf_tKp_y',{Lest_},{Gf_tKp_y_});
% get_Tf_u_ = casadi.Function('get_Tf_u',{Lest_},{Tf_u_});

% make anonymous functions
% get_Gf_tKp_u  = @(Lest) full(get_Gf_tKp_u_(Lest));
% get_Gf_tKp_y  = @(Lest) full(get_Gf_tKp_y_(Lest));
% get_Tf_u  = @(Lest) full(get_Tf_u_(Lest));

par = [up_(:); yp_(:); Lest_(:); yrf_(:); urf_(:)];

% construct cost function
er_y = yrf_(:) - yf_(:);
er_u = urf_(:) - uf_(:);
duf = uf_ - [up_(:,end) uf_(:,1:end-1)];                         % u_k - u_{k-1}
cost = er_y.'*Q*er_y + er_u(:).'*R*er_u(:) + duf(:).'*dR*duf(:); % total cost

% construct QP:
H = hessian(cost,uf_(:));
% get_H = casadi.Function('get_H',{par},{H});
cost_lin = cost - 0.5*uf_(:).'*H*uf_(:);
ct = jacobian(cost_lin,uf_(:));
get_ct = casadi.Function('get_ct',{par,uf_},{ct});
ct = get_ct(par,zeros(nu,f));
% get_ct = casadi.Function('get_ct',{par},{ct});
uf_sol = -H\ct.';
u_sol = uf_sol(1:nu,1);

% get solver
calc_u_v1 = casadi.Function('get_uk_analytic',{par},{u_sol});
calc_u_v2 = @(up,yp,Lest,yrf,urf) full(calc_u_v1([up(:); yp(:); Lest(:); yrf(:); urf(:)]));

% break down solver into contributions
% -> contribution from up
up2uk_ = jacobian(u_sol,up_(:));
up2uk_v1 = casadi.Function('get_upfac',{Lest_},{up2uk_});
up2usol = @(Lest) full(up2uk_v1(Lest));

% -> contribution from yp
yp2uk_ = jacobian(u_sol,yp_(:));
yp2uk_v1 = casadi.Function('get_ypfac',{Lest_},{yp2uk_});
yp2usol = @(Lest) full(yp2uk_v1(Lest));

% -> contribution from yrf
yrf2uk_ = jacobian(u_sol,yrf_(:));
yrf2uk_v1 = casadi.Function('get_yrffac',{Lest_},{yrf2uk_});
yrf2usol = @(Lest) full(yrf2uk_v1(Lest));

% -> contribution from urf
urf2uk_ = jacobian(u_sol,urf_(:));
urf2uk_v1 = casadi.Function('get_urffac',{Lest_},{urf2uk_});
urf2usol = @(Lest) full(urf2uk_v1(Lest));

% put results into a structure
get_usol.up  = up2usol;
get_usol.yp  = yp2usol;
get_usol.urf = urf2usol;
get_usol.yrf = yrf2usol;

if nargout > 1
% get Uf from Up, Uf, Yp, Yf, YRf URf
ICu_ = casadi.SX.zeros(nu*f,nu*(p+f));
ICy_ = casadi.SX.zeros(nu*f,ny*(p+f));
ICyr_ = casadi.SX.zeros(nu*f,ny*(2*f-1)); % YRf -> Uf
ICur_ = casadi.SX.zeros(nu*f,nu*(2*f-1)); % URf -> Uf
for k = 1:f
    ICu_((k-1)*nu+1:k*nu,(k-1)*nu+1:(k+p-1)*nu) = up2uk_;
    ICy_((k-1)*nu+1:k*nu,(k-1)*ny+1:(k+p-1)*ny) = yp2uk_;
    ICyr_((k-1)*nu+1:k*nu,(k-1)*ny+1:(k+f-1)*ny) = yrf2uk_;
    ICur_((k-1)*nu+1:k*nu,(k-1)*nu+1:(k+f-1)*nu) = urf2uk_;
end
ICup_ = ICu_(:,1:nu*p);     % Up -> Uf
ICuf_ = ICu_(:,nu*p+1:end); % Uf -> Uf
ICyp_ = ICy_(:,1:ny*p);     % Up -> Uf
ICyf_ = ICy_(:,ny*p+1:end); % Yf -> Uf

% Uf2Uf0_ = ICuf_ + ICyf_*Tf_u_;
% left_Uf = eye(nu*f) - Uf2Uf0_;
% Up2Uf_ = left_Uf\(ICup_+ICyf_*Gf_tKp_u_);
% Yp2Uf_ = left_Uf\(ICyp_+ICyf_*Gf_tKp_y_);
% YRf2Uf_ = left_Uf\ICyr_;
% URf2Uf_ = left_Uf\ICur_;
% inv_Uf = left_Uf\[ICup_+ICyf_*Gf_tKp_u_ ICyp_+ICyf_*Gf_tKp_y_ ICr_];
% Up2Uf_ = inv_Uf(:,1:size(ICup_,2));
% Yp2Uf_ = inv_Uf(:,size(ICup_,2)+1:size(ICup_,2)+size(ICyp_,2));
% Rf2Uf_ = inv_Uf(:,end-size(ICr_,2)+1:end);

% make functions out of contributions to Uf
% Up2Uf_v1 = casadi.Function('get_Up2Uf',{Lest_},{Up2Uf_});
% get_Up2Uf = @(Lest) full(Up2Uf_v1(Lest));
% Yp2Uf_v1 = casadi.Function('get_Yp2Uf',{Lest_},{Yp2Uf_});
% get_Yp2Uf = @(Lest) full(Yp2Uf_v1(Lest));
% YRf2Uf_v1 = casadi.Function('get_YRf2Uf',{Lest_},{YRf2Uf_});
% get_YRf2Uf = @(Lest) full(YRf2Uf_v1(Lest));
% URf2Uf_v1 = casadi.Function('get_URf2Uf',{Lest_},{URf2Uf_});
% get_URf2Uf = @(Lest) full(URf2Uf_v1(Lest));

% make functions out of contributions to Uf
ICup_v1 = casadi.Function('get_ICup',{Lest_},{ICup_});
get_ICup = @(Lest) full(ICup_v1(Lest));
ICuf_v1 = casadi.Function('get_ICuf',{Lest_},{ICuf_});
get_ICuf = @(Lest) full(ICuf_v1(Lest));
ICyp_v1 = casadi.Function('get_ICyp',{Lest_},{ICyp_});
get_ICyp = @(Lest) full(ICyp_v1(Lest));
ICyf_v1 = casadi.Function('get_ICyf',{Lest_},{ICyf_});
get_ICyf = @(Lest) full(ICyf_v1(Lest));
ICyr_v1 = casadi.Function('get_ICyr',{Lest_},{ICyr_});
get_ICyr = @(Lest) full(ICyr_v1(Lest));
ICur_v1 = casadi.Function('get_ICur',{Lest_},{ICur_});
get_ICur = @(Lest) full(ICur_v1(Lest));

% define optional outputs
varargout{1} = get_ICuf;
varargout{2} = get_ICyf;
varargout{3} = get_ICup;
varargout{4} = get_ICyp;
varargout{5} = get_ICyr;
varargout{6} = get_ICur;
varargout{7} = calc_u_v2;
end
end