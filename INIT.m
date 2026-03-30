%% INIT - Initial setup for the IVopt-DDPC project
%  This script:
%   1. Adds src/ and subdirectories to the MATLAB path
%   2. Configures cluster settings for Slurm job submission (optional)
%   3. Downloads and installs CasADi v3.6.7 binaries (needed for running MC simulations)
%   4. Downloads Crameri scientific colour maps (needed for plotting)

%% Adding paths
[pdir, ~  , ~] = fileparts(which(mfilename)); % find project directory
cd(pdir);
addpath(genpath(fullfile(pdir,'src')));       % add src & subdirectories to path

%% Specify cluster settings
% if you use a (Slurm) cluster, enter account, partition, and profile name of the cluster below, otherwise
% leave these as empty strings
ProfileName = ''; % Enter profile name used for the Slurm cluster (e.g., 'MyClusterProfile')
account     = ''; % Enter account e.g. 'institution-department'
partition   = ''; % Enter partitions to use, comma-separated if multiple (e.g., 'partition1,partition2')

% save these settings in src\SlurmSettings.mat
save(fullfile('src','SlurmSettings.mat'),'account','partition','ProfileName');

%% Downloading CasADi v3.6.7 binaries
% =========================================================================
%  Select CasADi v3.6.7 binaries & downloads them if not already present
%  Minimum required: R2018b (all platforms), R2023b (Mac Apple Silicon)
% ==========================================================================

CASADI_VERSION = '3.6.7';
destDir = fullfile(pwd, 'bin', sprintf('casadi-v%s', CASADI_VERSION));

if exist(destDir, 'dir') && ~isempty(dir(fullfile(destDir, '**')))
    fprintf('CasADi already installed, skipping download:\n  %s\n', destDir);
else
    CASADI_TAG     = '3.6.7';          % GitHub release tag
    BASE_URL       = sprintf('https://github.com/casadi/casadi/releases/download/%s/', CASADI_TAG);

    % --- 1. Detect MATLAB release --------------------------------------------
    matlabRelease = version('-release');          % e.g. '2023b'     %2024b
    matlabYear    = str2double(matlabRelease(1:end-1));
    matlabLetter  = lower(matlabRelease(end));      % 'a' or 'b'

    % Convert release to a comparable scalar: year + 0.5 if letter is 'b'
    releaseNum = matlabYear + (0.5 * strcmp(matlabLetter, 'b'));
    MIN_RELEASE        = 2018.5;   % R2018b
    MIN_RELEASE_APPLE  = 2023.5;   % R2023b (Apple Silicon native)

    fprintf('MATLAB release : %s\n', matlabRelease);

    % --- 2. Detect OS and architecture ---------------------------------------
    if ispc
        platform = 'windows';

    elseif ismac
        platform = 'mac';
        [status, arch] = system('uname -m');
        if status ~= 0
            error('Failed to query Mac architecture via ''uname -m''.');
        end
        macArch = strtrim(arch);   % 'arm64' = Apple Silicon, 'x86_64' = Intel
        fprintf('Mac architecture: %s\n', macArch);

    elseif isunix
        platform = 'linux';

    else
        error('Unrecognised operating system.');
    end

    fprintf('Platform       : %s\n', platform);

    % --- 3. Select URL -------------------------------------------------------
    switch platform

        case 'windows'
            if releaseNum < MIN_RELEASE
                error('CasADi %s for MATLAB on Windows requires R2018b or later. Detected: %s.', ...
                      CASADI_VERSION, matlabRelease);
            end
            filename = sprintf('casadi-%s-windows64-matlab2018b.zip', CASADI_VERSION);

        case 'linux'
            if releaseNum < MIN_RELEASE
                error('CasADi %s for MATLAB on Linux requires R2018b or later. Detected: %s.', ...
                      CASADI_VERSION, matlabRelease);
            end
            filename = sprintf('casadi-%s-linux64-matlab2018b.zip', CASADI_VERSION);

        case 'mac'
            if strcmp(macArch, 'arm64')
                % Apple Silicon: prefer native binary (R2023b+),
                % fall back to Rosetta binary (R2018b+)
                if releaseNum >= MIN_RELEASE_APPLE
                    filename = sprintf('casadi-%s-osx_arm64-matlab2018b.zip', CASADI_VERSION);
                    fprintf('Mac Apple Silicon: using native R2023b binary.\n');
                elseif releaseNum >= MIN_RELEASE
                    filename = sprintf('casadi-%s-osx64-matlab2018b.zip', CASADI_VERSION);
                    fprintf('Mac Apple Silicon: R2023b not met, falling back to Rosetta binary.\n');
                else
                    error('CasADi %s for MATLAB on Mac requires R2018b or later. Detected: %s.', ...
                          CASADI_VERSION, matlabRelease);
                end
            else
                % Intel Mac (x86_64)
                if releaseNum < MIN_RELEASE
                    error('CasADi %s for MATLAB on Mac (Intel) requires R2018b or later. Detected: %s.', ...
                          CASADI_VERSION, matlabRelease);
                end
                filename = sprintf('casadi-%s-osx64-matlab2018b.zip', CASADI_VERSION);
                fprintf('Mac Intel: using classic binary.\n');
            end
    end

    casadiUrl = [BASE_URL filename];
    fprintf('Download URL   : %s\n', casadiUrl);

    % --- 4. Download and unzip -----------------------------------------------
    destFile = fullfile(pwd, 'bin', filename);
    downloadAndUnzip(casadiUrl, destFile, destDir, 'CasADi');
end

% --- 5. Add to path and verify -------------------------------------------
addpath(destDir);
fprintf('CasADi added to MATLAB path.\n');

% Quick sanity check
try
    x = casadi.SX.sym('x');
    fprintf('CasADi %s loaded successfully.\n', CASADI_VERSION);
catch
    warning('CasADi was downloaded but could not be initialised. Check the path.');
end

%% Download Crameri colour binaries
% Crameri, Fabio. Scientific Colour Maps. Zenodo, 2019, doi:10.5281/ZENODO.1243862.

% --- 1. Download and unzip binaries -------------------------------------------
crameriUrl = 'https://nl.mathworks.com/matlabcentral/mlc-downloads/downloads/988c00c0-6131-42e7-8246-e1efbf8825a2/1da1cc21-d678-4dfc-9ab7-667916cc48ed/packages/zip';
destFile = fullfile(pwd, 'bin', 'crameri_colours.zip');
destDir  = fullfile(pwd, 'bin', 'crameri_colours');
downloadAndUnzip(crameriUrl, destFile, destDir, 'Crameri colour binaries');

% --- 2. Add to path and verify -------------------------------------------
addpath(destDir);
fprintf('Crameri colour binaries added to MATLAB path.\n');

try
    RGB = crameri('buda',5); % test if Crameri function is available
    fprintf('Crameri colour binaries loaded successfully.\n');
catch
    warning('Crameri colour binaries were downloaded but could not be initialised. Check the path.');
end

%% Local functions
function downloadAndUnzip(url, zipFile, extractDir, name)
    % Download, unzip, and clean up a zip file
    % Skips download if extractDir already exists.

    if exist(extractDir, 'dir') && ~isempty(dir(fullfile(extractDir, '**')))
        fprintf('Nonempty %s already exists, skipping download:\n  %s\n', name, extractDir);
        return
    end

    fprintf('Downloading %s ...\n', name);
    try
        websave(zipFile, url);
    catch ME
        error('Download failed: %s\nURL: %s', ME.message, url);
    end

    fprintf('Extracting to  : %s\n', extractDir);
    unzip(zipFile, extractDir);
    delete(zipFile);
    fprintf('Done.\n');
end