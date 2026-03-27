function get_usol = get_solver(opts,Q,R,dR)
%% Derives analytical feedback control laws for finite-horizon quadratic predictive control
%
% PROBLEM FORMULATION:
%   Solves the finite-horizon, unconstrained optimal control problem:
%       min_{u_f} ||y_{r,f} - y_f||²_Q + ||u_{r,f} - u_f||²_R + ||?u_f||²_dR
%   where:
%     - y_f = future outputs predictions based on Lf, past data, and future inputs
%     - y_{r,f}, u_{r,f} = reference trajectories
%     - ?u_f = change in future inputs 
%
% SOLUTION METHOD:
%   Derives the analytical solution using CasADi symbolic computations:
%       u_k = Cup(Lf)·up + Cyp(Lf)·yp + Curf(Lf)·urf + Cyrf(Lf)·yrf
%   where Cup, Cyp, Curf, Cyrf are matrix functions that are derived here.
%
% OUTPUT STRUCTURE (get_usol):
%   A struct with handles to functions computing control coefficients:
%     .up(Lf):   Contribution matrix Cup (feedback from past inputs)
%     .yp(Lf):   Contribution matrix Cyp (feedback from past outputs)
%     .urf(Lf):  Contribution matrix Curf (feedforward from future input reference)
%     .yrf(Lf):  Contribution matrix Cyrf (feedforward from future output reference)

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

% parameter vector
par = [up_(:); yp_(:); Lest_(:); yrf_(:); urf_(:)];

% construct cost function
er_y = yrf_(:) - yf_(:);
er_u = urf_(:) - uf_(:);
duf = uf_ - [up_(:,end) uf_(:,1:end-1)];                         % u_k - u_{k-1}
cost = er_y.'*Q*er_y + er_u(:).'*R*er_u(:) + duf(:).'*dR*duf(:); % total cost

% construct QP:
H = hessian(cost,uf_(:));
cost_lin = cost - 0.5*uf_(:).'*H*uf_(:);
ct = jacobian(cost_lin,uf_(:));
get_ct = casadi.Function('get_ct',{par,uf_},{ct});
ct = get_ct(par,zeros(nu,f));
uf_sol = -H\ct.';
u_sol = uf_sol(1:nu,1);

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

end