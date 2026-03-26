function m3_mat = process_m3_seed(cost_tot,cost_u,cost_y,Cases,cost_types)
%PROCESS_M3_SEED Compute per-seed DDPC performance costs for all cases
%   Returns a matrix numel(cost_types) x numel(Cases)

num_cost_types = numel(cost_types);
num_Cases = numel(Cases);
m3_mat = zeros(num_cost_types, num_Cases);

for kType = 1:num_cost_types
    cost_type = cost_types{kType};
    switch cost_type
        case 'cost_u'
            cost = cost_u;
        case 'cost_y'
            cost = cost_y;
        otherwise
            cost = cost_tot;
    end

    for kIVn = 1:num_Cases
        IVn = Cases{kIVn};
        m3_mat(kType, kIVn) = cost.(IVn);
    end
end
end
