function m1 = m1_data2struct(m1_UYf, Uf_ivs, pctiles)
num_Uf_ivs = numel(Uf_ivs);

fprintf("Processing m1 data\n")
m1 = struct();
tic
% Uf part
for k = 1:num_Uf_ivs
    data_k = squeeze( m1_UYf(k,:,:) );
    m1.Uf.(Uf_ivs{k}).data     = data_k;
    m1.Uf.(Uf_ivs{k}).mean     = mean(    data_k, 2);
    m1.Uf.(Uf_ivs{k}).median   = median(  data_k, 2);
    m1.Uf.(Uf_ivs{k}).pctiles  = prctile( data_k, pctiles, 2);
end
toc

end