function [SimCases,CaseDescr,noSimCases] = findCases2Sim(SimCases,noSimCases)
% Determine which IV method cases to simulate based on user specifications
% Handles flexible case selection with support for 'all', inclusion lists, and exclusion lists.
%
% INPUT:
%   SimCases    - cell array or cell containing 'all' string. If 'all', includes all available cases.
%                 Otherwise, only specified cases are included.
%   noSimCases  - cell array of cases to explicitly exclude from simulation.
%                 Applied AFTER inclusion logic (works with both specific lists and 'all' option).
%
% OUTPUT:
%   SimCases    - final cell array of cases to simulate (sorted)
%   CaseDescr   - descriptions of the cases to simulate
%   noSimCases  - final cell array of cases to exclude (complements SimCases)
%
% LOGIC:
%   1. If 'all' in SimCases -> start with all available cases
%   2. Otherwise -> use only cases specified in SimCases
%   3. Remove any cases in noSimCases
%   4. Return remaining cases as SimCases, complement set as noSimCases
%
% EXAMPLES:
%   [cases, noSim] = findCases2Sim({'all'}, {'iv4b','iv4c'});  % all except iv4b, iv4c
%   [cases, noSim] = findCases2Sim({'iv1','iv2','CLSPC'}, {});  % only iv1, iv2, CLSPC
%   [cases, noSim] = findCases2Sim({'iv1'}, {'iv1'});           % ERROR: no cases to simulate

% Get all available cases from CaseDefinitions
allCases = CaseDefinitions({});

% Step 1: Determine inclusion set based on SimCases specification
if ~isempty(intersect(SimCases,{'all'}))
    % 'all' keyword found -> start with all available cases
    SimCases = allCases;
    % noSimCases can still exclude specific cases from 'all'
end

% Step 2: Validate that all specified cases are valid (exist in CaseDefinitions)
checkCases(SimCases);
checkCases(noSimCases);

% Step 3: Compute final exclusion set = (cases not in SimCases) UNION (explicitly excluded)
noSimCases = union(setdiff(allCases,SimCases),noSimCases);

% Step 4: Compute final simulation set = SimCases minus exclusions
[SimCases,CaseDescr,~] = CaseDefinitions(noSimCases);

% Step 5: Error handling - ensure at least one case remains
if isempty(SimCases)
    error('No cases to simulate: the combination of Cases and noSimCases results in an empty set')
end
end

function checkCases(Cases)
    % Validate that all case names are recognized and defined in CaseDefinitions
    allCases = CaseDefinitions({});
    unknown = setdiff(Cases,allCases);
    if ~isempty(unknown)
        % Format unknown cases for error message
        unknownStr = sprintf(', %s', unknown{:});
        unknownStr(1:2) = [];  % remove leading ', '
        error('Unknown case(s) encountered: %s. Valid cases are: %s', unknownStr, sprintf('%s ', allCases{:}));
    end
end