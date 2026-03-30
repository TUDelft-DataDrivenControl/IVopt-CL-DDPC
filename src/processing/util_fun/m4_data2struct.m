function m4 = m4_data2struct(Yf_RelErr,Cases,nX)

fprintf("Processing m4 data\n")
m4 = struct;
tic;
for kC = 1:numel(Cases)
    caseName = Cases{kC};
    for iX = 1:nX
        iXstr = sprintf('iX%d',iX);
        % (Yf_hat  - Yf_hatS)./Yf_hatS -> yfhat
        m4.yfhat.(caseName).(iXstr).mean = squeeze(Yf_RelErr(iX,:,kC,2)); % mean over spX seeds
        m4.yfhat.(caseName).(iXstr).std  = squeeze(Yf_RelErr(iX,:,kC,1)); % mean std. dev. per seed over spX seeds
        m4.yfhat.(caseName).(iXstr).rms  = squeeze(Yf_RelErr(iX,:,kC,3)); % mean rms per seed over spX seeds
    end
end
toc

end