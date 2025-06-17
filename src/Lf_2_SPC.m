function [Cz_SPC,x0] = Lf_2_SPC(Lf,up2usol,yp2usol,urf2usol,yrf2usol,nu,ny,p,f,options)
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

nargoutchk(1,2);
if nargout == 2 && any(structfun(@isempty, options))
    error('Specify up, yp, urf, yrf to determine initial state of controller');
end

% ---------- create state-space representation of SPC controller ----------
% controller: u_k = Cu*up + Cy*yp + Cur*urf + Cyr*yrf
Cu  = up2usol(Lf);
Cy  = yp2usol(Lf);
Cur = urf2usol(Lf);
Cyr = yrf2usol(Lf);

C = [Cu Cy Cur Cyr];
A = blkdiag(eye(nu*(p-1)),zeros(nu),... up
            eye(ny*(p-1)),zeros(ny),... yp
            eye(nu*(f-1)),zeros(nu),... urf
            eye(ny*(f-1)),zeros(ny)); % yrf
A(nu*(p-1)+1:p*nu,:) = C;
B = [zeros(p*nu,2*ny+nu);...                    up
     blkdiag([zeros(ny*(p-1),ny);eye(ny)],...   yp
     [zeros(nu*(f-1),nu);eye(nu)],...           urf
     [zeros(nu*(f-1),ny);eye(nu)])];          % yrf
D = zeros(nu,2*ny+nu);

Cz_SPC = ss(A,B,C,D,[]);

% ----------------------- create initial state ----------------------------
if nargout == 2
    x0 = [options.up(:); ...
          options.yp(:); ...
          options.urf(:);...
          options.yrf(:)];
end

end