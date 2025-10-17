function m4 = m4_data2struct(Yf_RelErr,Cases,nX)

fprintf("Processing m4 data\n")
m4 = struct;
tic;
for kC = 1:numel(Cases)
    caseName = Cases{kC};
    for iX = 1:nX
        iXstr = sprintf('iX%d',iX);
        % (Yf_hat  - Yf)./Yf -> yfhat
        m4.yfhat.(caseName).(iXstr).mean = squeeze(Yf_RelErr(iX,:,kC,1,2)); % mean over spX seeds
        m4.yfhat.(caseName).(iXstr).std  = squeeze(Yf_RelErr(iX,:,kC,1,1)); % mean std. dev. per seed over spX seeds

        % (Yf_hatS  - Yf)./Yf -> yfhatS
        m4.yfhatS.(caseName).(iXstr).mean = squeeze(Yf_RelErr(iX,:,kC,2,2)); % mean over spX seeds
        m4.yfhatS.(caseName).(iXstr).std  = squeeze(Yf_RelErr(iX,:,kC,2,1)); % mean std. dev. per seed over spX seeds

        % Yf_by_Ep./Yf -> MatEp        ratio of contribution of Ep to Yf
        m4.MatEp.(caseName).(iXstr).mean = squeeze(Yf_RelErr(iX,:,kC,3,2)); % mean over spX seeds
        m4.MatEp.(caseName).(iXstr).std  = squeeze(Yf_RelErr(iX,:,kC,3,1)); % mean std. dev. per seed over spX seeds

        % Yf_by_Ef./Yf -> HfEf          = Hf*Ef./Yf
        m4.HfEf.(caseName).(iXstr).mean = squeeze(Yf_RelErr(iX,:,kC,4,2));  % mean over spX seeds
        m4.HfEf.(caseName).(iXstr).std  = squeeze(Yf_RelErr(iX,:,kC,4,1));  % mean std. dev. per seed over spX seeds
    end
end
toc

end