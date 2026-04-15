function [Cases,Descr] = CaseDefinitions(noSimCases)
%% Returns all simulated controller cases and their descriptions
% INPUT:
%   noSimCases: cell array of case names to exclude from the simulation (e.g. {'iv2b','iv3a'})

arguments
    noSimCases cell = {}
end

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
%   TrPred  | Transient Predictor

% make structure array for data
Cases = {'iv1','iv2a','iv2b','iv2c','iv3a','iv3c','iv4a','iv4b','iv4c',...
         'iv5a','iv5b','iv5c','iv6a','iv6c','CLSPC','actLf','TrPred'};
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
'SPC using the actual matrix Lf',...                      actLf  + SPC
'Transient Predictor'};%                                  TrPred

% make a selection of simulation cases
if ~isempty(noSimCases)
    [~,ia] = setdiff(Cases,noSimCases); ia = sort(ia);
    Cases = Cases(ia);
    Descr = Descr(ia);
end

end