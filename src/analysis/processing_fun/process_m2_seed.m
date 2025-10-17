function m2_mat = process_m2_seed(FroIDerror, Cases, IDerrorTypes)
%PROCESS_M2_SEED Compute identification error matrix for a single seed
%   m2_mat is size numel(IDerrorTypes) x numel(Cases)

num_IDerrorTypes = numel(IDerrorTypes);
num_Cases = numel(Cases);
m2_mat = zeros(num_IDerrorTypes, num_Cases);

for kType = 1:num_IDerrorTypes
    IDerrorType = IDerrorTypes{kType};
    for kIVn = 1:num_Cases
        IVn = Cases{kIVn};
        m2_mat(kType, kIVn) = FroIDerror.(IVn).(IDerrorType);
    end
end
end
