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
[Uf_2sls,Uf_iv4] = approx_IV_no_controller_info(u0,y0,yr0,opts);

% ---------------------- (5) w/ controller info. --------------------------
Uf_iv5 = approx_IV_controller_info(u0,y0,yr0,opts,Cz0);

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
        case 'iv2a'         % 2a) w/o Yf_iv (result for minimum asymptotic variance)
            Z0 = Uf_iv2;
    
% ---------------------- (3) LCF-IV ---------------------------------------
        case 'iv3a'         % 3a) orignal form of LCF-IV
            Z0 = [IV_Theta; Rf_yr0];
    
% ---------------------- (4) w/o controller info. -------------------------
        case 'iv4a'         % 4a) without Yf_iv
            Z0 = Uf_iv4;
    
% ---------------------- (5) w/ controller info. --------------------------
        case 'iv5a'         % 5a) without Yf_iv
            Z0 = Uf_iv5;
    
% ---------------------- (6) basic IV -------------------------------------
        case 'iv6a'         % 6a) IV composed of future reference signal
            Z0 = Rf_yr0;

% ---------------------- (7) IVopt + 2SLS ---------------------------------
        case 'iv7'         % 7) IV 2SLS applied to exogenous components
            Z0 = Uf_2sls;

% ---------------------- (8) iv4a but row by row --------------------------
        case 'iv8'         % 7) IV 2SLS applied to exogenous components
            Z0 = Uf_2sls;

% ---------------------- no IV to be made ---------------------------------
        otherwise
            continue;
    
    end
    
    % create IV based on Z0, TSLS flag, and description
    Z.add_IV(IV_name,Z0,descr=IV_descr);

end
end