clear; 
close all;
data_type = 'p'; % 'N', 'Re', or 'p' represented by X below
overwrite = false;

%% navigate to data\sys#\ref0_<>\dX\<subdir1>
[subdir1,src_dir] = get_subdir1(data_type);
cd(subdir1); % move to subdir1
fprintf('Analyzing data in %s\n',subdir1(length(src_dir)-3:end));

%% load data for choice of dX trials (as specified by <subdir1>)

% load dX settings from data\sys#\ref0_<>\dX\<subdir1>\dX_settings.mat
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
if isfile('processed_data.mat') && ~overwrite
    expectedVars = {'m1','mLf','m4'};
    availVars = {whos('-file','processed_data.mat').name};
    missingVars = setdiff(expectedVars, availVars);
    
    if isempty(missingVars)
        % do noting, all expected processed data is available
        fprintf('All expected processed data is available\n');

    else
        % add missing processed data
        fprintf('Not all necessary processed data found, adding missing data to file\n')
        process_dX(data_type,seeds,X_all,opts,plant,missingVars); % available variables will be emptyd
    end

else
    if ~isfile('processed_data.mat')
        % file not found -> create it
        fprintf("Processed data file 'processed_data.mat' not found\n");
        fprintf('Processing data in directory\n');
    else
        % file found, but needs to be overwritten
        fprintf('Overwriting processed_data.mat\n');
    end
    process_dX(data_type,seeds,X_all,opts,plant);
end

%% remove data path again
cd(src_dir);
rmpath(genpath(subdir1));