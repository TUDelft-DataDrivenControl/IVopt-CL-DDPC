function [Cases,Descr,Uf_ivs] = CaseDefinitions(noSimCases)
%% Returns all simulated controller cases and their descriptions
% INPUT:
%   noSimCases: cell array of case names to exclude from the simulation (e.g. {'iv2','iv3'})

arguments
    noSimCases cell = {}
end

%    Name   |       Description
% ----------|----------------------------------------------------------
%   iv1     | open-loop IV
%   iv2     | opt IV
%   iv3     | LCF-IV
%   iv4a    | opt IV w/o Cz0 info: basic (2SLS)
%   iv4b    | opt IV w/o Cz0 info: 2SLS + causal + time invariant
%   iv4c    | opt IV w/o Cz0 info: row by row estimation
%   iv4d    | opt IV w/o Cz0 info: row by row + time invariant
%   iv5b    | opt IV w/o Cz0 info: basic (2SLS)
%   iv6     | basic IV: future reference
%   CLSPC   | CL-SPC
%   actLf   | SPC using the actual matrix Lf
%   TrPred  | Transient Predictor
% methods with IVs use SPC-type of predictor

% make structure array for data
Cases = {'iv1','iv2','iv3','iv4a','iv4b','iv4c','iv4d',...
        'iv5a','iv5b','iv5c','iv5d','iv6','CLSPC','actLf','TrPred'};
Uf_ivs = setdiff(Cases,{'CLSPC','actLf','TrPred','iv3'}); % IVs akin to Uf
Descr = {...
'open-loop IV',...                                                          iv1    + SPC
'optimal IV',...                                                            iv2    + SPC
'LCF-IV',...                                                                iv3    + SPC
'approx. opt. IV w/o controller info: 2SLS',...                             iv4a   + SPC
'approx. opt. IV w/o controller info: 2SLS + lower-blk-diag-avged',...      iv4b   + SPC
'approx. opt. IV w/o controller info: row-wise est.',...                    iv4c   + SPC
'approx. opt. IV w/o controller info: row-wise est. + blk-diag-avged',...   iv4d   + SPC
'approx. opt. IV w/  controller info: 2SLS',...                             iv5a   + SPC
'approx. opt. IV w/  controller info: 2SLS + lower-blk-diag-avged',...      iv5b   + SPC
'approx. opt. IV w/  controller info: row-wise est.',...                    iv5c   + SPC
'approx. opt. IV w/  controller info: row-wise est. + blk-diag-avged',...   iv5d   + SPC
'basic IV: future reference',...                                            iv6    + SPC
'CL-SPC',...                                                                CLSPC
'SPC using the actual matrix Lf',...                                        actLf  + SPC
'Transient Predictor'};%                                                    TrPred

% make a selection of simulation cases
if ~isempty(noSimCases)
    [~,ia] = setdiff(Cases,noSimCases); ia = sort(ia);
    Cases = Cases(ia);
    Descr = Descr(ia);
    Uf_ivs = intersect(Uf_ivs,Cases);
end

end