function m0 = m0_data2struct(m0_UYf, Uf_ivs, Yf_ivs, nX)
num_Uf_ivs = numel(Uf_ivs);
num_Yf_ivs = numel(Yf_ivs);

fprintf("Processing m0 data\n")
m0 = struct;
% Uf part
tic
for k = 1:num_Uf_ivs
    for iX = 1:nX
        iXstr = sprintf('iX%d',iX);
        m0.Uf.(Uf_ivs{k}).(iXstr).mean     = m0_UYf{k,iX,1};
        m0.Uf.(Uf_ivs{k}).(iXstr).median   = m0_UYf{k,iX,2};
        m0.Uf.(Uf_ivs{k}).(iXstr).pctiles  = m0_UYf{k,iX,3};
    end
end
toc

% Yf part
tic
for k = num_Uf_ivs+(1:num_Yf_ivs)
    kY = k - num_Uf_ivs;
    for iX = 1:nX
        iXstr = sprintf('iX%d',iX);
        m0.Yf.(Yf_ivs{kY}).(iXstr).mean     = m0_UYf{k,iX,1};
        m0.Yf.(Yf_ivs{kY}).(iXstr).median   = m0_UYf{k,iX,2};
        m0.Yf.(Yf_ivs{kY}).(iXstr).pctiles  = m0_UYf{k,iX,3};
    end
end
toc

end