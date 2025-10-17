function m4 = m4_data2struct(Yf_RelErr_mean,Yf_RelErr_sd,Cases,nX)

fprintf("Processing m4 data\n")
m4 = struct;
tic;
for kC = 1:numel(Cases)
    caseName = Cases{kC};
    for iX = 1:nX
        iXstr = sprintf('iX%d',iX);
        % (Yf_hat  - Yf)./Yf -> yfhat
        m4.yfhat.(caseName).(iXstr).mean = squeeze(Yf_RelErr_mean(iX,:,kC,1)); % mean over spX seeds
        m4.yfhat.(caseName).(iXstr).std  = squeeze(Yf_RelErr_sd(iX,:,kC,1));   % mean std. dev. per seed over spX seeds

        % (Yf_hatS  - Yf)./Yf -> yfhatS
        m4.yfhatS.(caseName).(iXstr).mean = squeeze(Yf_RelErr_mean(iX,:,kC,2)); % mean over spX seeds
        m4.yfhatS.(caseName).(iXstr).std  = squeeze(Yf_RelErr_sd(iX,:,kC,2));   % mean std. dev. per seed over spX seeds

        % Yf_by_Ep./Yf -> MatEp        ratio of contribution of Ep to Yf
        m4.MatEp.(caseName).(iXstr).mean = squeeze(Yf_RelErr_mean(iX,:,kC,3)); % mean over spX seeds
        m4.MatEp.(caseName).(iXstr).std  = squeeze(Yf_RelErr_sd(iX,:,kC,3));   % mean std. dev. per seed over spX seeds

        % Yf_by_Ef./Yf -> HfEf          = Hf*Ef./Yf
        m4.HfEf.(caseName).(iXstr).mean = squeeze(Yf_RelErr_mean(iX,:,kC,4)); % mean over spX seeds
        m4.HfEf.(caseName).(iXstr).std  = squeeze(Yf_RelErr_sd(iX,:,kC,4));   % mean std. dev. per seed over spX seeds
    end
end
toc

end