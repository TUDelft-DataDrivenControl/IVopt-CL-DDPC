function [Cz_SPC, x0] = Lf_2_SPC(Lf, up2usol, yp2usol, urf2usol, yrf2usol, nu, ny, p, f, options)
arguments (Input)
    Lf      double
    up2usol function_handle
    yp2usol function_handle
    urf2usol function_handle
    yrf2usol function_handle
    nu       (1,1) double
    ny       (1,1) double
    p        (1,1) double
    f        (1,1) double
    options.up double = []
    options.yp double = []
    options.urf double = []
    options.yrf double = []
end
arguments (Output)
    Cz_SPC ss
    x0 double
end

% --- Get controller matrices ---
Cup  = up2usol(Lf);
Cyp  = yp2usol(Lf);
Curf = urf2usol(Lf);
Cyrf = yrf2usol(Lf);
C = [Cup Cyp Curf Cyrf];  % u_k = C * x_k

% --- Sizes ---
n_up   = p * nu;
n_yp   = p * ny;
n_urf  = f * nu;
n_yrf  = f * ny;
nx = n_up + n_yp + n_urf + n_yrf;
n_input = ny + nu + ny; % [y_k; ur_{k+f}; yr_{k+f}]

% --- Construct shift matrices ---
% Shift up (u_{k-p} to u_{k-1}), insert u_k at the bottom
A_up = [zeros((p-1)*nu, nu), eye((p-1)*nu)];
A_up = [A_up; zeros(nu, n_up)];

% Shift yp (y_{k-p} to y_{k-1}), insert y_k at bottom via B
A_yp = [zeros((p-1)*ny, ny), eye((p-1)*ny)];
A_yp = [A_yp; zeros(ny, n_yp)];

% Shift urf (ur_k to ur_{k+f-1}), insert ur_{k+f} at bottom via B
A_urf = [zeros((f-1)*nu, nu), eye((f-1)*nu)];
A_urf = [A_urf; zeros(nu, n_urf)];

% Shift yrf (yr_k to yr_{k+f-1}), insert yr_{k+f} at bottom via B
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

% --- Initial state (optional) ---
if nargout == 2
    if any(structfun(@isempty, options))
        error('Specify up, yp, urf, yrf to determine initial state of controller');
    end
    x0 = [options.up(:); options.yp(:); options.urf(:); options.yrf(:)];
end
end
