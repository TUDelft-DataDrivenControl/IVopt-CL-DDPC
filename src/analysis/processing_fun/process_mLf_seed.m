function mLf_data_iX_ks = process_mLf_seed(Cases,ny,f,cols_Lf,Lf,ks)
%PROCESS_MLF_SEED Compute per-seed identified Lf for all cases
%   Returns a 3D matrix of size num_Cases x (ny*f) x cols_Lf
num_Cases = numel(Cases);

mLf_data_iX_ks = zeros(num_Cases,ny*f,cols_Lf);
for kC = 1:num_Cases
    CaseName = Cases{kC};
    if ~(strcmp(CaseName,'actLf') && ks > 1)
        mLf_data_iX_ks(kC,:,:) = Lf.(CaseName);
    end
end

end