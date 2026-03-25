clear; 
close all;
data_type = 'Re'; % 'N', 'Re', or 'p' represented by X below
overwrite = false;
plotting = false;

if ismember('SlurmProfile1',parallel.clusterProfiles) % running on cluster?
    onCluster = true;
else
    onCluster = false;
end

%% navigate to data\raw\sys#\ref0_<>\dX\<subdir1>
[subdir1,src_dir] = get_subdir1(data_type);
cd(subdir1); % move to subdir1
fprintf('Analyzing data in %s\n',subdir1(length(src_dir)-3:end));

%% load data for choice of dX trials (as specified by <subdir1>)

% load dX settings from data\raw\sys#\ref0_<>\dX\<subdir1>\dX_settings.mat
load(sprintf('d%s_settings.mat',data_type));

switch data_type
    case 'N'
        spX = spN;
        [p,f,nu,ny] = deal(opts.p,opts.f,opts.nu,opts.ny);
        X_all = N_all;
    case 'Re'
        spX = spRe;
        [p,f,nu,ny,N] = deal(opts.p,opts.f,opts.nu,opts.ny,opts.N);
        X_all = Re_all;
    case 'p'
        spX = spP;
        [f,nu,ny,N] = deal(opts.f,opts.nu,opts.ny,opts.N);
        X_all = p_all;
end
nX = numel(X_all);

% get data structures - load or by processing data
if isfile('processed_data.mat')
    expectedVars = {'m0','m1','m2','mLf','m3','m4'};
    availVars = {whos('-file','processed_data.mat').name};
    missingVars = setdiff(expectedVars, availVars);

    if isempty(missingVars)
        % all processed data available, load it
        fprintf('All necessary processed data found, loading file\n')
        load("processed_data.mat");

    elseif ~isempty(missingVars) && overwrite
        % overwrite existing processed data
        fprintf('Not all necessary processed data found, overwriting file\n')
        [m0,m1,m2,mLf,m3,m4] = process_dX(data_type,seeds,X_all,opts,plant);
        
    else
        % add only missing processed data
        fprintf('Not all necessary processed data found, adding missing data to file\n')
        [m0,m1,m2,mLf,m3,m4] = process_dX(data_type,seeds,X_all,opts,plant,missingVars); % available variables will be empty
        clear(availVars{:}); % clear empty available variables before loading them
        load("processed_data.mat",availVars{:}); % load existing variables
    end

else
    fprintf('Processed data file not found\n')
    fprintf('Processing data in directory\n');
    [m0,m1,m2,mLf,m3,m4] = process_dX(data_type,seeds,X_all,opts,plant);
end

%% remove data path again
cd(src_dir);
rmpath(genpath(subdir1));