function [Z,nlcf]= get_Z(u0,y0,yr0,Cz0,Uf_iv2,Cases,Descr,opts)
%% creates all of the IVs using the IV_4_DDPC class
% The main objective of this function is to create an instance of the
% IV_4_DDPC class in which to save the created IVs

[p,f,DMCS] = deal(opts.p,opts.f,opts.DMCS);
nCz = numel(Cases);

%% ========================= calculate IVs ================================
% some necessary calculations before assigning IVs using IV_4_DDPC class
fprintf('Obtaining IVs...\n');

% create data matrices
Uf_r01 = make_Page(u0(:,p+1:end),f,DMCS); % - data w/ noise
Yf_r01 = make_Page(y0(:,p+1:end),f,DMCS); % - data w/ noise

% ---------------------- (3) LCF-IV ---------------------------------------
% see "Data-Driven Predictive Control Using  Closed-Loop Data: An
% Instrumental  Variable Approach" (2023) by Wang et al.
% DOI: 10.1109/LCSYS.2023.3340444
[~,Vc,Uc] = lncf(ss(Cz0));
Hcv = make_blk_tril_toeplitz(Vc.A,Vc.B,Vc.C,Vc.D,f);
Hcu = make_blk_tril_toeplitz(Uc.A,Uc.B,Uc.C,Uc.D,f);
IV_Theta = Hcv*Uf_r01 + Hcu*Yf_r01;
nlcf = size(Vc.C,1); % needed to get IV_Theta from Z later
% -> also, if nlcf = ny, lets interpret IV_Theta as an IV resembling noiseless future outputs

Rf_yr0 = make_Page(yr0(:,p+1:end),f,DMCS); % also used for (6)

% ---------------------- (4) w/o controller info. -------------------------
% w/o -> don't known Cz0 exactly, but know rho & feedback configuration
[Uf_iv4a,Uf_iv4b,Uf_iv4c,Uf_iv4d] = approx_IV_no_controller_info(u0,y0,yr0,opts,Cases);

% ---------------------- (5) w/ controller info. --------------------------
[Uf_iv5a,Uf_iv5b,Uf_iv5c,Uf_iv5d] = approx_IV_controller_info(u0,y0,yr0,opts,Cz0,Cases);

%% ========================== assign IVs ==================================
% using an instance of the IV_4_DDPC class

% ---------------------- (1) open-loop IV ---------------------------------
% 1) i.e. 'normal' least-squares regression)
Z = IV_4_DDPC(u0,y0,p,f,DMCS); % initializes IV object & makes open-loop IV ('iv1')

for kIV = 1:nCz % loop over Cases (which contains IV names)
    IV_name = Cases{kIV};    % IV name
    IV_descr = Descr{kIV};   % IV description

    switch IV_name
% ---------------------- (2) optimal IV -----------------------------------
        case 'iv2'         % 2) w/o Yf_iv (result for minimum asymptotic variance)
            Z0 = Uf_iv2;
    
% ---------------------- (3) LCF-IV ---------------------------------------
        case 'iv3'         % 3) orignal form of LCF-IV
            Z0 = [IV_Theta; Rf_yr0];      
    
% ---------------------- (4) w/o controller info. ----------------------
        case 'iv4a'        % 4a) 2SLS applied to exogenous components
            Z0 = Uf_iv4a;

        case 'iv4b'        % 4b) 2SLS + causality + time invariance (averaged blk-diags)
            Z0 = Uf_iv4b;
    
        case 'iv4c'        % 4c) row-by-row approach = causal by construction
            Z0 = Uf_iv4c;
    
        case 'iv4d'        % 4d) row-by-row + time invariance
            Z0 = Uf_iv4d;
    
% ---------------------- (5) w/ controller info. --------------------------
        case 'iv5a'        % 5a) 2SLS applied to exogenous components
            Z0 = Uf_iv5a;

        case 'iv5b'        % 5b) 2SLS + causality + time invariance (averaged blk-diags)
            Z0 = Uf_iv5b;
    
        case 'iv5c'        % 5c) row-by-row approach = causal by construction
            Z0 = Uf_iv5c;
    
        case 'iv5d'        % 5d) row-by-row + time invariance
            Z0 = Uf_iv5d;
    
% ---------------------- (6) basic IV -------------------------------------
        case 'iv6'         % 6a) IV composed of future reference signal
            Z0 = Rf_yr0;

% ---------------------- no IV to be made ---------------------------------
        otherwise
            continue;
    
    end
    
    % create IV based on Z0, TSLS flag, and description
    Z.add_IV(IV_name,Z0,descr=IV_descr);

end
end