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
%   iv3a    | SPC using LCF-IV
%   iv4a    | SPC using approx. opt. IV w/o controller info.
%   iv5a    | SPC using approx. opt. IV w/  controller info.
%   iv6a    | SPC using basic IV: future reference
%   iv7     | SPC using optimal IV components + 2SLS
%   CLSPC   | CL-SPC
%   actLf   | SPC using the actual matrix Lf
%   TrPred  | Transient Predictor

% make structure array for data
Cases = {'iv1','iv2a','iv3a','iv4a','iv5a','iv6a','iv7','CLSPC','actLf','TrPred'};
Descr = {...
'open-loop IV',...                                        iv1    + SPC
'optimal IV',...                                          iv2a   + SPC
'LCF-IV',...                                              iv3a   + SPC
'approx. opt. IV w/o controller info.',...                iv4a   + SPC
'approx. opt. IV w/  controller info.',...                iv5a   + SPC
'basic IV: future reference',...                          iv6a   + SPC
'IVopt: 2SLS',...                                         iv7    + SPC
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