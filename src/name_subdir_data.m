function subdir = name_subdir_data(opts)
[Re, p, f, N, Ncl, dRk, Rk, Qk] = deal(opts.Re, opts.p, opts.f, opts.N, opts.Ncl, opts.dRk, opts.Rk, opts.Qk);

% Helper function to trim to minimal digits in scientific notation
trimmed_exp = @(x) regexprep(sprintf('%e', x), '(\.\d*?)0+(e[+-]?\d+)', '$1$2'); % trims trailing 0s
% Also remove . if nothing follows
trimmed_exp = @(x) regexprep(trimmed_exp(x), '\.(e)', '$1');

% Apply formatting
Re_str  = trimmed_exp(Re); N_str   = trimmed_exp(N);  Ncl_str = trimmed_exp(Ncl);
Qk_str  = trimmed_exp(Qk); Rk_str  = trimmed_exp(Rk); dRk_str = trimmed_exp(dRk);
subdir = sprintf('Re_%s_p_%d_f_%d_N_%s_Ncl_%s_Qk_%s_Rk_%s_dRk_%s',Re_str,p,f,N_str,Ncl_str,Qk_str,Rk_str,dRk_str);
subdir = replace(subdir,'.','p');
subdir = replace(subdir,'+','');
end