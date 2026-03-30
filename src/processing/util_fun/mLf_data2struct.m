function mLf = mLf_data2struct(mLf_data, Cases, nX, pctiles)

fprintf("Processing mLf data\n")
mLf = struct(); % <- mLf_data = nan(nX,num_Cases,spX,ny*f,size(Lf,1),size(Lf,2));
tic
for kC = 1:numel(Cases)
    caseName = Cases{kC};
    for iX = 1:nX
        iXstr = sprintf('iX%d',iX);
        mLf_data_iX_kC = squeeze(mLf_data{iX}(kC, :, :,:));
        if strcmp(caseName,'actLf')
            mLf.(caseName).(iXstr) = squeeze(mLf_data_iX_kC(1,:,:)); % same for all seeds
        else
            mLf.(caseName).(iXstr).data    = mLf_data_iX_kC;
            mLf.(caseName).(iXstr).mean    = squeeze( mean(    mLf_data_iX_kC, 1) );
            mLf.(caseName).(iXstr).median  = squeeze( median(  mLf_data_iX_kC, 1) );
            mLf.(caseName).(iXstr).pctiles = squeeze( prctile( mLf_data_iX_kC, pctiles, 1) );
            mLf.(caseName).(iXstr).std     = squeeze( std( mLf_data_iX_kC, 0 ,1) );
        end
    end
end
toc

end