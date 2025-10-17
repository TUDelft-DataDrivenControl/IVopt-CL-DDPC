function m3 = m3_data2struct(m3_data,cost_types,Cases,pctiles)

fprintf("Processing m3 data\n")
m3 = struct(); % <- m3_data(kType, kIVn, iX, ks)
tic
for kType = 1:numel(cost_types)
    costName = cost_types{kType};
    for kIVn = 1:numel(Cases)
        caseName = Cases{kIVn};
        m3.(costName).(caseName).data    = squeeze( m3_data(kType, kIVn, :, :) );
        m3.(costName).(caseName).mean    = mean(    m3.(costName).(caseName).data, 2);
        m3.(costName).(caseName).median  = median(  m3.(costName).(caseName).data, 2);
        m3.(costName).(caseName).pctiles = prctile( m3.(costName).(caseName).data, pctiles, 2);
    end
end
toc

end