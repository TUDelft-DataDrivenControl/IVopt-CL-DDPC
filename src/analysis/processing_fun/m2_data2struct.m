function m2 = m2_data2struct(m2_data, Cases, IDerrorTypes, pctiles)

fprintf("Processing m2 data\n")
m2 = struct(); % <- m2_data(kType, kIVn, iX, ks)
tic
for kType = 1:numel(IDerrorTypes)
    typeName = IDerrorTypes{kType};
    for kIVn = 1:numel(Cases)
        caseName = Cases{kIVn};
        m2.(typeName).(caseName).data    = squeeze( m2_data(kType, kIVn, :, :) );
        m2.(typeName).(caseName).mean    = mean(    m2.(typeName).(caseName).data, 2);
        m2.(typeName).(caseName).median  = median(  m2.(typeName).(caseName).data, 2);
        m2.(typeName).(caseName).pctiles = prctile( m2.(typeName).(caseName).data, pctiles, 2);
    end
end
toc

end