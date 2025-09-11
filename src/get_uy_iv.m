function [u_iv,y_iv] = get_uy_iv(Z,nlcf,yr0,opts)
[nu,ny,p,f] = deal(opts.nu,opts.ny,opts.p,opts.f);
iyf = nu*f+(1:ny*f);    % indices to select Yf from Z

% ---------------------- (1) open-loop IV ---------------------------------
% see u0 & y0, std. dev. = 0

% ---------------------- (2) optimal IV -----------------------------------
% 2a) w/o Yf_iv
Uf_iv2 = Z.iv2a_;
[u_iv.m2a,u_iv.s2a] = diag_stats(Uf_iv2,nr=nu,anti=true);

% 2b) w/  Yf_iv
Yf_iv2 = Z.iv2b_(iyf,:);
[y_iv.m2b,y_iv.s2b] = diag_stats(Yf_iv2,nr=ny,anti=true);

% 3a) w/  Yf_iv + 2SLS
Uf_iv2c = Z.iv2c_;                                        
[u_iv.m2c,u_iv.s2c] = diag_stats(Uf_iv2c,nr=nu,anti=true); 

% ---------------------- (3) LCF-IV ---------------------------------------
% 2a)          -> interpreting IV_Theta as IV for outputs
IV_Theta = Z.iv3a_(1:nlcf*f,:);
[y_iv.m3a,y_iv.s3a] = diag_stats(IV_Theta,nr=nlcf,anti=true);

% 2b) w/ 2SLS
Uf_iv3c = Z.iv3c_;
[u_iv.m3c,u_iv.s3c] = diag_stats(Uf_iv3c,nr=nu,anti=true);

% ---------------------- (4) w/o controller info. -------------------------
% 4a) w/o Yf_iv
Uf_iv4 = Z.iv4a_;
[u_iv.m4a,u_iv.s4a] = diag_stats(Uf_iv4,nr=nu,anti=true);

% 4b) w/  Yf_iv
Yf_iv4 = Z.iv4b_(iyf,:);
[y_iv.m4b,y_iv.s4b] = diag_stats(Yf_iv4,nr=ny,anti=true);

% 4c) w/  Yf_iv + 2SLS
Uf_iv4c = Z.iv4c_;
[u_iv.m4c,u_iv.s4c] = diag_stats(Uf_iv4c,nr=nu,anti=true);

% ---------------------- (5) w/ controller info. --------------------------
% 5a) w/o Yf_iv
Uf_iv5 = Z.iv5a_;
[u_iv.m5a,u_iv.s5a] = diag_stats(Uf_iv5,nr=nu,anti=true);

% 5b) w/  Yf_iv
Yf_iv5 = Z.iv5b_(iyf,:);
[u_iv.m5b,u_iv.s5b] = diag_stats(Yf_iv5,nr=ny,anti=true);

% 5c) w/  Yf_iv + 2SLS
Uf_iv5c = Z.iv5c_;
[u_iv.m5c,u_iv.s5c] = diag_stats(Uf_iv5c,nr=nu,anti=true);

% ---------------------- (6) basic IV --------------------------
% 6a) reference
y_iv.m6a = yr0(:,p+1:end); 

% 6c) reference + 2SLS
Uf_iv6c = Z.iv6c_;
[u_iv.m6c,u_iv.s6c] = diag_stats(Uf_iv6c,nr=nu,anti=true);
end