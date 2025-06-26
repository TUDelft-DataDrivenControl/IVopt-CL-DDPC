fig10 = figure(10); 
Fig10FS = 12;
tl10 = tiledlayout(4+num_iter*1,4,'TileSpacing','compact','Padding','compact');

% Row 1  - plot real system parameters
nexttile;
nexttile;
nexttile;
imagesc_vik([Up2Yf_inno Uf2Yf_inno]);
nexttile;
imagesc_vik(Yp2Yf_inno);

% Row 2
nexttile;
imagesc_parts(Ef_r1, uRp, uRf, Yp_r1);
title('$\frac{1}{N} E_f \mathcal{Z}^\top$', 'Interpreter', 'latex', 'FontSize', Fig10FS);
ylabel('$\mathcal{Z}=[R_p; R_f; Y_p]$', 'Interpreter', 'latex', 'FontSize', Fig10FS);
nexttile;
imagesc_parts(Wp_r1, uRp, uRf, Yp_r1);
title('$\frac{1}{N} W_p \mathcal{Z}^\top$, $W_p=[U_p; U_f; Y_p]$', 'Interpreter', 'latex', 'FontSize', Fig10FS);
nexttile; 
LestZr = Yf_r1 * Ziv.' / N * pinv(Wp_r1 * Ziv.' / N);
imagesc_vik(LestZr(:,1:end-ny*p));
title('$\frac{1}{N} Y_f\mathcal{Z}^\top (\frac{1}{N} W_p\mathcal{Z}^\top)^\dagger$, $U_{pf} \rightarrow Y_f$', 'Interpreter', 'latex', 'FontSize', Fig10FS);
nexttile;
imagesc_vik(LestZr(:,end-ny*p+1:end));
title('$\frac{1}{N} Y_f\mathcal{Z}^\top (\frac{1}{N} W_p\mathcal{Z}^\top)^\dagger$, $Y_p \rightarrow Y_f$', 'Interpreter', 'latex', 'FontSize', Fig10FS);

% Row 3
nexttile;
imagesc_parts(Ef_r1, Up_r1, Uf_r1, Yp_r1);
ylabel('$\mathcal{Z}=[U_p; U_f; Y_p]$', 'Interpreter', 'latex', 'FontSize', Fig10FS);
nexttile;
imagesc_parts(Wp_r1, Up_r1, Uf_r1, Yp_r1);
nexttile; 
LestWp = Yf_r1 * Wp_r1.' / N * pinv(Wp_r1 * Wp_r1.' / N);
imagesc_vik(LestWp(:,1:end-ny*p));
nexttile;
imagesc_vik(LestWp(:,end-ny*p+1:end));

% Row 4
nexttile;
imagesc_parts(Ef_r1, Up_r2, Uf_r2, Yp_r2);
ylabel('$\mathcal{Z}=[\tilde{U}_p; \tilde{U}_f; \tilde{Y}_p]_{\mathrm{opt}}$', 'Interpreter', 'latex', 'FontSize', Fig10FS);
nexttile;
imagesc_parts(Wp_r1, Up_r2, Uf_r2, Yp_r2);
nexttile; 
LestZopt = Yf_r1 * Ziv_opt.' / N * pinv(Wp_r1 * Ziv_opt.' / N);
imagesc_vik(LestZopt(:,1:end-ny*p));
nexttile;
imagesc_vik(LestZopt(:,end-ny*p+1:end));

%%
% produce imagesc of correlation
function imagesc_parts(Ef,varargin)
narginchk(4,inf);
Z = [];
nvarargin = length(varargin);
xvals = zeros(1,nvarargin);
for k = 1:nvarargin
    Z = [Z;varargin{k}];
    xvals(k) = size(varargin{k},1);
end
corEfZ = Ef*Z.'/size(Ef,2);
imagesc_vik(corEfZ);
xvals = cumsum(xvals)+0.5;
for k = 1:nvarargin-1
    xline(xvals(k),'LineWidth',1.5);
end
end