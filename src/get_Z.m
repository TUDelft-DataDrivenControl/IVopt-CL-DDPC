function [Z,Cases,nlcf]= get_Z(u0,y0,yr0,Cz0,Uf_iv2,Yf_iv2,opts)
%% creates all of the IVs using the IV_4_DDPC class
% The main objective of this function is to create an instance of the
% IV_4_DDPC class in which to save the created IVs

[p,f] = deal(opts.p,opts.f);

%% ============= create data structure and explain cases ==================
%    Name   |       Description
% ----------|----------------------------------------------------------
%   iv1     | SPC using open-loop IV
%   iv2a    | SPC using optimal IV
%   iv2b    | SPC using optimal IV + Yf_iv
%   iv2c    | SPC using optimal IV + Yf_iv + 2SLS
%   iv3a    | SPC using LCF-IV
%   iv3c    | SPC using LCF-IV + 2SLS
%   iv4a    | SPC using approx. opt. IV w/o controller info.
%   iv4b    | SPC using approx. opt. IV w/o controller info. + Yf_iv
%   iv4c    | SPC using approx. opt. IV w/o controller info. + Yf_iv + 2SLS
%   iv5a    | SPC using approx. opt. IV w/  controller info.
%   iv5b    | SPC using approx. opt. IV w/  controller info. + Yf_iv
%   iv5c    | SPC using approx. opt. IV w/  controller info. + Yf_iv + 2SLS
%   iv6a    | SPC using basic IV: future reference
%   iv6c    | SPC using basic IV: future reference + 2SLS
%   CLSPC   | CL-SPC
%   actLf   | SPC using the actual matrix Lf

% make structure array for data
Cases = {'iv1','iv2a','iv2b','iv2c','iv3a','iv3c','iv4a','iv4b','iv4c',...
         'iv5a','iv5b','iv5c','iv6a','iv6c','CLSPC','actLf'};
Descr = {...
'open-loop IV',...                                        iv1    + SPC
'optimal IV',...                                          iv2a   + SPC
'optimal IV + Yf_iv',...                                  iv2b   + SPC
'optimal IV + Yf_iv + 2SLS',...                           iv2c   + SPC
'LCF-IV',...                                              iv3a   + SPC
'LCF-IV + 2SLS',...                                       iv3c   + SPC
'approx. opt. IV w/o controller info.',...                iv4a   + SPC
'approx. opt. IV w/o controller info. + Yf_iv',...        iv4b   + SPC
'approx. opt. IV w/o controller info. + Yf_iv + 2SLS',... iv4c   + SPC
'approx. opt. IV w/  controller info.',...                iv5a   + SPC
'approx. opt. IV w/  controller info. + Yf_iv',...        iv5b   + SPC
'approx. opt. IV w/  controller info. + Yf_iv + 2SLS',... iv5c   + SPC
'basic IV: future reference',...                          iv6a   + SPC   
'basic IV: future reference + 2SLS',...                   iv6c   + SPC
'CL-SPC',...                                              CLSPC
'SPC using the actual matrix Lf'};%                       actLf  + SPC
nCz = numel(Cases);

%% ========================= calculate IVs ================================
% some necessary calculations before assigning IVs using IV_4_DDPC class
fprintf('Obtaining IVs...\n');

% create data matrices
[~,~,Uf_r01] = make_Hankel(u0,p,f); % - data w/ noise
[~,~,Yf_r01] = make_Hankel(y0,p,f);

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

[~,~,Rf_yr0] = make_Hankel(yr0,p,f); % also used for (6)

% ---------------------- (4) w/o controller info. -------------------------
% w/o -> don't known Cz0 exactly, but know rho & feedback configuration
[Uf_iv4,Yf_iv4] = approx_IV_no_controller_info(u0,y0,yr0,opts);

% ---------------------- (5) w/ controller info. --------------------------
[Uf_iv5,Yf_iv5] = approx_IV_controller_info(u0,y0,yr0,opts,Cz0);

%% ========================== assign IVs ==================================
% using an instance of the IV_4_DDPC class

% ---------------------- (1) open-loop IV ---------------------------------
% 1) i.e. 'normal' least-squares regression)
Z = IV_4_DDPC(u0,y0,p,f); % initializes IV object & makes open-loop IV ('iv1')

for kIV = 2:nCz-2 % loop over remaining IV names
    IV_name = Cases{kIV};    % IV name
    IV_descr = Descr{kIV};   % IV description

    switch IV_name
% ---------------------- (2) optimal IV -----------------------------------
        case 'iv2a'         % 2a) w/o Yf_iv (result for minimum asymptotic variance)
            Z0 = Uf_iv2;
            method_flag = 0;
        case 'iv2b'         % 2b) = iv2a + Yf_iv
            Z0 = [Uf_iv2;Yf_iv2];
            method_flag = 0;
        case 'iv2c'         % 2c) = iv2a + Yf_iv + 2SLS
            Z0 = [Uf_iv2;Yf_iv2];
            method_flag = 1;
    
% ---------------------- (3) LCF-IV ---------------------------------------
        case 'iv3a'         % 3a) orignal form of LCF-IV
            Z0 = [IV_Theta; Rf_yr0];
            method_flag = 0;
        case 'iv3c'         % 3c) => 3a) + 2SLS
            Z0 = [IV_Theta; Rf_yr0];
            method_flag = 1;
    
% ---------------------- (4) w/o controller info. -------------------------
        case 'iv4a'         % 4a) without Yf_iv
            Z0 = Uf_iv4;
            method_flag = 0;
        case 'iv4b'         % 4b) => 4a) + Yf_iv
            Z0 = [Uf_iv4;Yf_iv4];
            method_flag = 0;
        case 'iv4c'         % 4c) => 4a) + Yf_iv + 2SLS
            Z0 = [Uf_iv4;Yf_iv4];
            method_flag = 1;
    
% ---------------------- (5) w/ controller info. --------------------------
        case 'iv5a'         % 5a) without Yf_iv
            Z0 = Uf_iv5;
            method_flag = 0;
        case 'iv5b'         % 5b) => 5a) + Yf_iv
            Z0 = [Uf_iv5;Yf_iv5];
            method_flag = 0;
        case 'iv5c'         % 5c) => 5a) + Yf_iv + 2SLS
            Z0 = [Uf_iv5;Yf_iv5];
            method_flag = 1;
    
% ---------------------- (6) basic IV -------------------------------------
        case 'iv6a'         % 6a) IV composed of future reference signal
            Z0 = Rf_yr0;
            method_flag = 0;
        case 'iv6c'         % 6c) => 6a) + SLS
            Z0 = Rf_yr0;
            method_flag = 1;
    
    end
    
    % create IV based on Z0, TSLS flag, and description
    Z.add_IV(IV_name,Z0,method=method_flag,descr=IV_descr);

end
end