function [Cz_SPC, x0] = Lf_2_SPC(Lf, usol_funs, opts, sigs, options)
%% Constructs a state-space SPC controller from learned dynamics
%
% STATE-SPACE STRUCTURE:
%   x_k = [u_{k-p:k-1}; y_{k-p:k-1}; u_{r,k:k+f-1}; y_{r,k:k+f-1}]
%   
%   - Segment 1 (up): p past inputs u_{k-p} to u_{k-1}
%   - Segment 2 (yp): p past outputs y_{k-p} to y_{k-1}
%   - Segment 3 (urf): f future reference inputs u_r (u_{r,k:k+f-1})
%   - Segment 4 (yrf): f future reference outputs y_r (y_{r,k:k+f-1})
%
%   A: Block-diagonal shift operators that progressively update the past/future history.
%      Each block shifts data and inserts new values. The last row of the up segment
%      implements the control law: u_k = C*x_k (computed from Lf via usol_funs).
%
%   B: Maps 3 inputs to update the state:
%      - y_k:        updates past outputs segment (yp)
%      - u_{r,k+f}:  updates future reference inputs (urf)
%      - y_{r,k+f}:  updates future reference outputs (yrf)
%
%   C: Extracts control input u_k from state using the data-driven control law.
%      C is computed as [Cup, Cyp, Curf, Cyrf] where each term encodes the
%      dependency of u_k on past inputs, past outputs, future reference inputs,
%      and future reference outputs respectively.
%
%   D: Zero matrix (no feedthrough).
%
% INPUT ARGUMENTS:
%   Lf:        An estimated matrix mapping p past inputs/outputs and f future inputs
%              to f future outputs. Size: (ny*f) × ((nu+ny)*p + nu*f)
%   usol_funs: Struct with functions mapping Lf to control law coefficients:
%              .up(Lf):  effect of past inputs on u_k
%              .yp(Lf):  effect of past outputs on u_k
%              .urf(Lf): effect of future reference inputs on u_k
%              .yrf(Lf): effect of future reference outputs on u_k
%   opts:      Struct containing nu (inputs), ny (outputs), p (past horizon), f (future horizon)
%   sigs:      Struct with signal names for labeling inputs/outputs
%   options:   Optional initial condition [up, yp, urf, yrf]

arguments (Input)
    Lf      double
    usol_funs struct
    opts     struct
    sigs struct % signal names
    options.up double = []
    options.yp double = []
    options.urf double = []
    options.yrf double = []
end
arguments (Output)
    Cz_SPC ss
    x0 double
end
[nu,ny,p,f] = deal(opts.nu,opts.ny,opts.p,opts.f);

% --- Get controller matrices ---
Cup  = usol_funs.up(Lf);   % up  -> u_k
Cyp  = usol_funs.yp(Lf);   % yp  -> u_k
Curf = usol_funs.urf(Lf);  % urf -> u_k
Cyrf = usol_funs.yrf(Lf);  % yrf -> u_k
C = [Cup Cyp Curf Cyrf];   % u_k = C * x_k

% --- Sizes ---
n_up   = p * nu;
n_yp   = p * ny;
n_urf  = f * nu;
n_yrf  = f * ny;
nx = n_up + n_yp + n_urf + n_yrf;
n_input = ny + nu + ny; % [y_k; u_{r,k+f}; y_{r,k+f}]

% --- Construct shift matrices ---
% Shift up (u_{k-p} to u_{k-1}), insert u_k at the bottom
A_up = [zeros((p-1)*nu, nu), eye((p-1)*nu)];
A_up = [A_up; zeros(nu, n_up)];

% Shift yp (y_{k-p} to y_{k-1}), insert y_k at bottom via B
A_yp = [zeros((p-1)*ny, ny), eye((p-1)*ny)];
A_yp = [A_yp; zeros(ny, n_yp)];

% Shift urf (ur_k to ur_{k+f-1}), insert u_{r,k+f} at bottom via B
A_urf = [zeros((f-1)*nu, nu), eye((f-1)*nu)];
A_urf = [A_urf; zeros(nu, n_urf)];

% Shift yrf (yr_k to yr_{k+f-1}), insert y_{r,k+f} at bottom via B
A_yrf = [zeros((f-1)*ny, ny), eye((f-1)*ny)];
A_yrf = [A_yrf; zeros(ny, n_yrf)];

% Block diagonal A
A = blkdiag(A_up, A_yp, A_urf, A_yrf);

% --- Insert u_k = C x_k into last row of up segment ---
A(n_up - nu + (1:nu), :) = C;

% --- Input mapping matrix B ---
B = zeros(nx, n_input);

% y_k updates last block of yp
B(n_up + n_yp - ny + (1:ny), 1:ny) = eye(ny);

% ur_{k+f} updates last block of urf
B(n_up + n_yp + n_urf - nu + (1:nu), ny + (1:nu)) = eye(nu);

% yr_{k+f} updates last block of yrf
B(n_up + n_yp + n_urf + n_yrf - ny + (1:ny), ny + nu + (1:ny)) = eye(ny);

% --- Output: u_k = C x_k ---
D = zeros(nu, n_input);

% --- Create state-space system ---
Cz_SPC = ss(A, B, C, D, []);

% naming the state-space system
Cz_SPC.u(1:ny)          = sigs.yk;   % y_k
Cz_SPC.u(ny+1:ny+nu)    = sigs.urkf; % u_{r,k+f}
Cz_SPC.u(end-ny+(1:ny)) = sigs.yrkf; % y_{r,k+f}
Cz_SPC.y                = sigs.uk;   % u_k

% --- Initial state (optional) ---
if nargout == 2
    if any(structfun(@isempty, options))
        error('Specify up, yp, urf, yrf to determine initial state of controller');
    end
    x0 = [options.up(:); options.yp(:); options.urf(:); options.yrf(:)];
end
end
