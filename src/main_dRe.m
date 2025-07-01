function main_dRe(nRe,iRe,opts,opts2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           DDPC using an Optimal-IV
%           Authors: R. Dinkla, T. Oomen, J.W. van Wingerden
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this is a wrapper function for main.m that provides a systematic way to
% iterate over Re
arguments (Input)
    % number of Re values in range [Re_min Re_max]
    nRe       (1,1) double {mustBeInteger,mustBeGreaterThanOrEqual(nRe,1)}
    % index of Re value in above range
    iRe       (1,1) double {mustBeInteger,mustBeGreaterThanOrEqual(iRe,1),mustBeLessThanOrEqual(iRe,nRe)} 
    opts.Re_min (1,1) double = 4.81e-5;            % min. innovation noise variance
    opts.Re_max (1,1) double = 4.81e-1;            % max. innovation noise variance
    opts2.plot       logical = false;
    opts2.p    (1,1) double  = 20;       % window lengths
    opts2.f    (1,1) double  = 20;
    opts2.N    (1,1) double  = 1e4;
    opts2.Ncl  (1,1) double  = 1500;     % simulation length of SPC
    opts2.dRk  (1,1) double  = 1;        % weights
    opts2.Rk   (1,1) double  = 1;
    opts2.Qk   (1,1) double  = 1e2;
    opts2.seed (1,1) double  = 1;
    opts2.save       logical = true;     % save data
end

%% get value of N
Re_all = logspace(log10(opts.Re_min),log10(opts.Re_max),nRe);
opts2.Re = Re_all(iRe);

%% ================= handle path in raw data directory ====================
% folder structure:
%   dRe \ Re_<Re_min>_<Re_max>_<nRe>_<the rest> \   Re_%s    \ seed_%d.mat
% = dRe \             <subdir1>                 \ <subdir2>  \ seed_%d.mat
opts2.raw_dir = {'dRe'};

% ------------------------- set name of subdir1 ---------------------------
[N, p, f, Ncl, dRk, Rk, Qk] = deal(opts2.N, opts2.p, opts2.f, opts2.Ncl, opts2.dRk, opts2.Rk, opts2.Qk);

% Helper function to trim to minimal digits in scientific notation
trimmed_exp = @(x) regexprep(sprintf('%e', x), '(\.\d*?)0+(e[+-]?\d+)', '$1$2'); % trims trailing 0s
% Also remove . if nothing follows
trimmed_exp = @(x) regexprep(trimmed_exp(x), '\.(e)', '$1');

% Apply formatting
Re_min_s = trimmed_exp(opts.Re_min);
Re_max_s = trimmed_exp(opts.Re_max);
N_s  = trimmed_exp(N); Ncl_s = trimmed_exp(Ncl);
Qk_s  = trimmed_exp(Qk); Rk_s  = trimmed_exp(Rk); dRk_s = trimmed_exp(dRk);
subdir1 = sprintf('Re_%s_%s_%d_p_%d_f_%d_N_%s_Ncl_%s_Qk_%s_Rk_%s_dRk_%s',...
                   Re_min_s,Re_max_s,nRe,p,f,N_s,Ncl_s, Qk_s, Rk_s, dRk_s);
subdir1 = replace(subdir1,'.','p');
subdir1 = replace(subdir1,'+','');
opts2.raw_dir{2} = subdir1;

% ------------------------- set name of subdir2 ---------------------------
nDigits = ceil(log10(nRe + 1));
format_subdir2 = ['%0', num2str(nDigits), 'd_Re_%.2e'];
subdir2 = sprintf(format_subdir2,iRe,opts2.Re);
subdir2 = replace(subdir2,'.','p'); % replace . with p
opts2.raw_dir{3} = subdir2;

%% ==================== pass options to main function ======================
opts2 = namedargs2cell(opts2);
main(opts2{:});
end