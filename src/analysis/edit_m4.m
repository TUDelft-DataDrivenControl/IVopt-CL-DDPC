clear; 
close all;
load("pdir.mat",'pdir'); % load path of project directory
src_dir = fullfile(pdir,'src');
data_dir = fullfile(pdir,'data');

% Search recursively for all 'process_data.mat' files
filePattern = fullfile(data_dir, '**', 'processed_data.mat');
fileList = dir(filePattern);

%% navigate to data\raw\sys#\ref0_<>\dX\<subdir1>
% iterate over folders with a 'processed_data.mat' file
for kFn = 1:length(fileList)
subdir1 = fileList(kFn).folder;
cd(subdir1); % move to subdir1

% check if edit has been made before
if isfile('log.mat')
   m4_edited = load("log.mat",'log').log.m4_edited;
   if m4_edited
       continue
   end
end

if contains(subdir1, ['dRe' filesep 'Re_'])
        data_type = 'Re';
elseif contains(subdir1, ['dN' filesep 'N_'])
        data_type = 'N';
elseif contains(subdir1, ['dp' filesep 'p_'])
        data_type = 'p';
else
    error('Data type not recognized');
end

fprintf('Rewriting m4 data in %s\n',subdir1(length(src_dir)-3:end));

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

if ~isfile('processed_data.mat')
    error('processed_data.mat file not found')
end

% recalculate m4 and save to processed_data.mat
[~,~,~,~,~,~] = process_dX(data_type,seeds,X_all,opts,plant,{'m4'});

% save \ update log.mat file
if isfile("log.mat")
    log = load("log.mat",'log').log;
end
log.m4_edited = true;
save("log.mat",'log');
end
cd(src_dir);