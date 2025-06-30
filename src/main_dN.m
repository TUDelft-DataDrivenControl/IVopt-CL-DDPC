function main_dN(nN,iN,opts,opts2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           DDPC using an Optimal-IV
%           Authors: R. Dinkla, T. Oomen, J.W. van Wingerden
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this is a wrapper function for main.m that provides a systematic way to
% iterate over N
arguments (Input)
    % number of N values in range [Nmin Nmax]
    nN        (1,1) double {mustBeInteger,mustBeGreaterThanOrEqual(nN,1)}
    % index of N value in above range
    iN        (1,1) double {mustBeInteger,mustBeGreaterThanOrEqual(iN,1),mustBeLessThanOrEqual(iN,nN)} 
    opts.Nmin (1,1) double = 100;            % min. number of Hankel matrix columns
    opts.Nmax (1,1) double = 1e5;            % max. number of Hankel matrix columns
    opts2.Re   (1,1) double  = 4.81e-2;  % innovation noise variance
    opts2.plot       logical = false;
    opts2.p    (1,1) double  = 20;       % window lengths
    opts2.f    (1,1) double  = 20;
    opts2.Ncl  (1,1) double  = 1500;     % simulation length of SPC
    opts2.dRk  (1,1) double  = 1;        % weights
    opts2.Rk   (1,1) double  = 1;
    opts2.Qk   (1,1) double  = 1e2;
    opts2.seed (1,1) double  = 1;
    opts2.save       logical = true;     % save data
end

%% get value of N
N_all = floor(logspace(log10(opts.Nmin),log10(opts.Nmax),nN));
opts2.N = N_all(iN);

%% ================= handle path in raw data directory ====================
% folder structure:
%   dN \ N_<Nmin>_<Nmax>_<nN>_<the rest> \   N_%d    \ seed_%d.mat
% = dN \             <subdir1>              \ <subdir2> \ seed_%d.mat
opts2.raw_dir = {'dN'};

% ------------------------- set name of subdir1 ---------------------------
[Re, p, f, Ncl, dRk, Rk, Qk] = deal(opts2.Re, opts2.p, opts2.f, opts2.Ncl, opts2.dRk, opts2.Rk, opts2.Qk);

% Helper function to trim to minimal digits in scientific notation
trimmed_exp = @(x) regexprep(sprintf('%e', x), '(\.\d*?)0+(e[+-]?\d+)', '$1$2'); % trims trailing 0s
% Also remove . if nothing follows
trimmed_exp = @(x) regexprep(trimmed_exp(x), '\.(e)', '$1');

% Apply formatting
Nmin_s = trimmed_exp(opts.Nmin);
Nmax_s = trimmed_exp(opts.Nmax);
Re_s  = trimmed_exp(Re); Ncl_s = trimmed_exp(Ncl);
Qk_s  = trimmed_exp(Qk); Rk_s  = trimmed_exp(Rk); dRk_s = trimmed_exp(dRk);
subdir1 = sprintf('N_%s_%s_%d_Re_%s_p_%d_f_%d_Ncl_%s_Qk_%s_Rk_%s_dRk_%s',...
                   Nmin_s,Nmax_s,nN,Re_s,p,f, Ncl_s, Qk_s, Rk_s, dRk_s);
subdir1 = replace(subdir1,'.','p');
subdir1 = replace(subdir1,'+','');
opts2.raw_dir{2} = subdir1;

% ------------------------- set name of subdir2 ---------------------------
nDigits = ceil(log10(nN + 1));
format_subdir2 = ['%0', num2str(nDigits), 'd_N_%d'];
subdir2 = sprintf(format_subdir2,iN,opts2.N);
opts2.raw_dir{3} = subdir2;

%% ==================== pass options to main function ======================
opts2 = namedargs2cell(opts2);
main(opts2{:});
end