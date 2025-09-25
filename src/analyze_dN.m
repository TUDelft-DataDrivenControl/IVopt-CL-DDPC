clear; 
close all;

%% navigate to data\raw\dN\<subdir1>
load("pdir.mat"); % load path of project directory
src_dir = fullfile(pdir,'src');
raw_dir = fullfile(pdir,'data','raw');
dN_dir = fullfile(pdir,'data','raw','dN');

% find <subdir1> candidates:
% -> subdirectories in data\raw\dN that match the naming convention from name_subdir1 in main_dN
subdir1s = find_named_subdirs(dN_dir);

% choose <subdir1> to use
subdir1 = choose_named_subdir(subdir1s);
subdir1 = fullfile(dN_dir,subdir1); % set path to subdir1
addpath(genpath(subdir1));

%% load data for choice of dN trials (as specified by <subdir1>)
cd(subdir1); % move to subdir1

% load dN settings from data\raw\dN\<subdir1>\dN_settings.mat
load('dN_settings.mat');
spX = spN;
[p,f,nu,ny] = deal(opts.p,opts.f,opts.nu,opts.ny);

if isfile('processed_data.mat')
    load("processed_data.mat");
else
    fprintf('Processed data file not found\n')
    fprintf('Processing data in directory\n');
    process_dN;
end

%% Plotting - figure 1: example of Uf_iv

pX = 5; % index of X to plot this for
opts.N = N_all(pX);

iX_str = sprintf('iX%d',pX);
% select the 25th and 75th percentile -> 6 & 16th of 21
% see size(m0.Uf.iv1.iX1.pctiles,3);

Uf_ivs = fieldnames(m0.Uf);
Uf_lb_fs = replace(Uf_ivs,'iv','l');
Uf_ub_fs = replace(Uf_ivs,'iv','u');
Uf_av_fs = replace(Uf_ivs,'iv','m');
num_Uf_ivs = numel(Uf_ivs);

u_iv = struct;
for kIV = 1:num_Uf_ivs
    lb_name = Uf_lb_fs{kIV};
    ub_name = Uf_ub_fs{kIV};
    av_name = Uf_av_fs{kIV};
    iv_name = Uf_ivs{kIV};
    
    u_iv.(av_name) = m0.Uf.(iv_name).(iX_str).median;
    u_iv.(lb_name) = squeeze(m0.Uf.(iv_name).(iX_str).pctiles(:,:,6));
    u_iv.(ub_name) = squeeze(m0.Uf.(iv_name).(iX_str).pctiles(:,:,16));
end

Yf_ivs = fieldnames(m0.Yf);
Yf_lb_fs = replace(Yf_ivs,'iv','l');
Yf_ub_fs = replace(Yf_ivs,'iv','u');
Yf_av_fs = replace(Yf_ivs,'iv','m');
num_Yf_ivs = numel(Yf_ivs);

y_iv = struct;
for kIV = 1:num_Yf_ivs
    lb_name = Yf_lb_fs{kIV};
    ub_name = Yf_ub_fs{kIV};
    av_name = Yf_av_fs{kIV};
    iv_name = Yf_ivs{kIV};
    
    y_iv.(av_name) = m0.Yf.(iv_name).(iX_str).median;
    y_iv.(lb_name) = squeeze(m0.Yf.(iv_name).(iX_str).pctiles(:,:,6));
    y_iv.(ub_name) = squeeze(m0.Yf.(iv_name).(iX_str).pctiles(:,:,16));
end
plot_IV_trajectories(u_iv,y_iv,opts,"useBounds",true);
clear u_iv y_iv

%% Plotting - figure 2: quality of optimal IV approx. vs. N
% -> difference of Uf_iv w.r.t. Uf_iv2a
% -> difference of Yf_iv w.r.t. Yf_iv2b

% plotting Uf frobenius norm errors
% select the 25th and 75th percentile -> 6 & 16th of 21
% see size(m1.Uf.iv1.pctiles,2);

% --------------------- set colours and make figure -----------------------
load("pdir.mat",'pdir');
addpath(fullfile(pdir,'bin','external','crameri_colours'))
cCram = crameri('roma', 9); % colors

fig2 = figure();
fig2.Units = 'pixels';
fig2.Position = [2600 500 1000 600];
tl2  = tiledlayout(2,1,"TileSpacing",'compact');
fontSize_xyLabel = 15;

% --------------------- plotting for Uf_ivs -------------------------------
ax2_1 = nexttile(tl2);

Uf_ivs     = fieldnames(m1.Uf);
Uf_lb_fs = replace(Uf_ivs,'iv','l');
Uf_ub_fs = replace(Uf_ivs,'iv','u');
Uf_av_fs = replace(Uf_ivs,'iv','m');
num_Uf_ivs = numel(Uf_ivs);

sqrt_Nall = diag(sqrt(N_all));
u_iv = struct;
for kIV = 1:num_Uf_ivs
    lb_name = Uf_lb_fs{kIV};
    ub_name = Uf_ub_fs{kIV};
    av_name = Uf_av_fs{kIV};
    iv_name = Uf_ivs{kIV};
    
    u_iv.(av_name) = sqrt_Nall\m1.Uf.(iv_name).median;
    u_iv.(lb_name) = sqrt_Nall\squeeze(m1.Uf.(iv_name).pctiles(:,6));
    u_iv.(ub_name) = sqrt_Nall\squeeze(m1.Uf.(iv_name).pctiles(:,16));
    
    % determine color
    switch iv_name
        case 'iv6c', label = '$j=6c$: ref + 2SLS';         col = [0 0 0];       LineStyle = '-';
        case 'iv2a', label = '$j=2a$: opt. IV';            col = cCram(1,:);    LineStyle = '-.';
        case 'iv2c', label = '$j=2c$: opt. IV + 2SLS';     col = cCram(2,:);    LineStyle = '-.';
        case 'iv3c', label = '$j=3c$: LCF + 2SLS';         col = cCram(3,:);    LineStyle = '--';
        case 'iv4a', label = '$j=4a$: w/o Cz info';        col = cCram(4,:);    LineStyle = ':';
        case 'iv4c', label = '$j=4c$: w/o Cz info + 2SLS'; col = cCram(5,:);    LineStyle = ':';
        case 'iv5a', label = '$j=5a$: w/ Cz info';         col = cCram(6,:);    LineStyle = '-.';
        case 'iv5c', label = '$j=5c$: w/ Cz info + 2SLS';  col = cCram(7,:);    LineStyle = '-.';
        case 'iv1',  label = '$j=1$: OL-IV';               col = cCram(9,:);    LineStyle = '-.';
        otherwise, label = iv_name;                    col = [0 0 0];       LineStyle = '-';
    end

    plotLineWithFill(u_iv.(av_name),u_iv.(lb_name),u_iv.(ub_name), label, ...
        Color=col, FaceAlpha = 0.25, LineStyle=LineStyle, LineWidth=2,...
         x=N_all, XScale='log');
    hold on;

    u_entries(kIV) = struct('Color',col, 'Alpha',0.25, ...
                    'LineStyle',LineStyle, 'LineWidth',2, ...
                    'Text',label);
end
grid on;
ylabel('$\frac{\|\Delta_j U_{\mathrm{f}}^{\mathrm{iv},2a}\|_2}{\sqrt{N}}$',...
    'Interpreter','latex','Rotation',0,'FontSize',fontSize_xyLabel*1.25);
axLeg_2_1 = customLegend(u_entries,ax2_1);

% --------------------- plotting for Yf_ivs -------------------------------
ax2_2 = nexttile(tl2);

Yf_ivs     = fieldnames(m1.Yf);
Yf_lb_fs = replace(Yf_ivs,'iv','l');
Yf_ub_fs = replace(Yf_ivs,'iv','u');
Yf_av_fs = replace(Yf_ivs,'iv','m');
num_Yf_ivs = numel(Yf_ivs);

y_iv = struct;
for kIV = 1:num_Yf_ivs
    lb_name = Yf_lb_fs{kIV};
    ub_name = Yf_ub_fs{kIV};
    av_name = Yf_av_fs{kIV};
    iv_name = Yf_ivs{kIV};
    
    y_iv.(av_name) = sqrt_Nall\m1.Yf.(iv_name).median;
    y_iv.(lb_name) = sqrt_Nall\squeeze(m1.Yf.(iv_name).pctiles(:,6));
    y_iv.(ub_name) = sqrt_Nall\squeeze(m1.Yf.(iv_name).pctiles(:,16));

    % determine color
    switch iv_name
        case 'iv1',  label = '$j=1$: OL-IV';         col = cCram(9,:);  LineStyle = '--';
        case 'iv2b', label = '$j=2b$: opt. IV';      col = cCram(1,:);  LineStyle = '-.';
        case 'iv4b', label = '$j=4b$: w/o Cz info';  col = cCram(4,:);  LineStyle = ':';
        case 'iv5b', label = '$j=5b$: w/ Cz info';   col = cCram(6,:);  LineStyle = '-.';
        case 'iv3a', label = '$j=3a$: LCF-IV Theta'; col = cCram(8,:);  LineStyle = '--';
        case 'iv6a', label = '$j=3a,6a$: ref';       col = [0 0 0];     LineStyle = '-';
        otherwise,   label = field;              col = [0 0 0];     LineStyle = '-';
    end

    plotLineWithFill(y_iv.(av_name),y_iv.(lb_name),y_iv.(ub_name), label, ...
        Color=col, FaceAlpha = 0.25, LineStyle=LineStyle, LineWidth=2,...
         x=N_all, XScale='log');
    hold on;

    y_entries(kIV) = struct('Color',col, 'Alpha',0.25, ...
                    'LineStyle',LineStyle, 'LineWidth',2, ...
                    'Text',label);
end
grid on;
ylabel('$\frac{\|\Delta_j Y_{\mathrm{f}}^{\mathrm{iv},2b}\|_2}{\sqrt{N}}$',...
    'Interpreter','latex','Rotation',0,'FontSize',fontSize_xyLabel*1.25)
xlabel('$N$','interpreter','latex','FontSize',fontSize_xyLabel);
axLeg_2_2 = customLegend(y_entries,ax2_2);

%% remove data path again
cd(src_dir);
rmpath(genpath(subdir1));

%% Helper functions
function subdirs = find_named_subdirs(parentDir)
% FIND_NAMED_SUBDIRS returns subdirectories in parentDir that match
% the naming convention from name_subdir1 in main_dN
%
% Example:
%   subdirs = find_named_subdirs(pwd);

    % List all directories in parentDir
    d = dir(parentDir);
    d = d([d.isdir]);             % keep only directories
    names = {d.name};
    
    % Remove '.' and '..'
    names = names(~ismember(names,{'.','..'}));
    
    % Regex pattern that matches the naming convention
    %   Example: N_1p0e0_1p0e1_50_sys_1_Re_1p0e3_p_2_f_3_Ncl_1p0e2_Qk_1p0e1_Rk_5p0e0_dRk_1p0e0
    pattern = ['^N_[-0-9ep]+_[-0-9ep]+_\d+_sys_\d+_' ...   % Nmin, Nmax, nN, sys
               'Re_[-0-9ep]+_p_\d+_f_\d+_' ...            % Re, p, f
               'Ncl_[-0-9ep]+_Qk_[-0-9ep]+_Rk_[-0-9ep]+_dRk_[-0-9ep]+$']; % Ncl, Qk, Rk, dRk
    
    % Keep only those that match
    isMatch = cellfun(@(x) ~isempty(regexp(x, pattern, 'once')), names);
    subdirs = names(isMatch);
end

function chosenDir = choose_named_subdir(subdirs)
% CHOOSE_NAMED_SUBDIR chooses a subdirectory from subdirs
    
if numel(subdirs) > 1
    while true
        % Display options
        fprintf('Available subdirectories:\n');
        for i = 1:numel(subdirs)
            fprintf('  [%d] %s\n', i, subdirs{i});
        end

        % Prompt user
        choice = input('Select a subdirectory by number: ', 's');

        % Validate choice
        choiceNum = str2double(choice);
        if ~isnan(choiceNum) && choiceNum >= 1 && choiceNum <= numel(subdirs)
            % Valid choice -> exit loop
            break;
        else
            clc;
            fprintf('Invalid choice. Please select a number between 1 and %d.\n\n', numel(subdirs));
        end
    end
    clc;
    % Get chosen directory (full path)
    chosenDir = subdirs{choiceNum};
else
    chosenDir = subdirs{1};
end

end

function h = plotLineWithFill(y, yLower, yUpper, description, varargin)
% plotLineWithFill plots a line with shaded bounds and a single legend
% entry that shows both the line and the patch together.
%
% Usage:
%   plotLineWithFill(y, yLower, yUpper, description)
%   plotLineWithFill(y, yLower, yUpper, description, 'x', xVals)
%   plotLineWithFill(y, yLower, yUpper, description, 'LineStyle','--','FaceAlpha',0.5)
%   plotLineWithFill(y, yLower, yUpper, description, 'XScale','log')
%
% Inputs:
%   y          - main line values
%   yLower     - lower bound (same length as y)
%   yUpper     - upper bound (same length as y)
%   description- string used for legend
%
% Optional name-value pairs:
%   'x'        - x values (default = 1:length(y))
%   'LineStyle' - line style for main line (default = '-')
%   'FaceAlpha' - transparency of fill (default = 0.25)
%   'Color'     - RGB triplet or color char (default = next color in ColorOrder)
%   'LineWidth' - line width (default = 1.5)
%   'XScale'    - 'linear' (default) or 'log'
%
% Output:
%   h - structure with handles to line and fill

    % Ensure column-to-row formatting
    y      = y(:).';
    yLower = yLower(:).';
    yUpper = yUpper(:).';

    % Parse inputs
    p = inputParser;
    addParameter(p, 'LineStyle', '-', @ischar);
    addParameter(p, 'FaceAlpha', 0.25, @(a)isnumeric(a)&&isscalar(a));
    addParameter(p, 'Color', [], @(c)(ischar(c) || (isnumeric(c)&&numel(c)==3)));
    addParameter(p, 'LineWidth', 1.5, @(a)isnumeric(a)&&isscalar(a));
    addParameter(p, 'x', 1:length(y), @isnumeric);
    addParameter(p, 'XScale', 'linear', @(s)ischar(s) && any(strcmpi(s,{'linear','log'})));
    parse(p, varargin{:});

    % Extract results
    ls     = p.Results.LineStyle;
    fa     = p.Results.FaceAlpha;
    c      = p.Results.Color;
    lw     = p.Results.LineWidth;
    x      = p.Results.x(:).';
    xScale = lower(p.Results.XScale);

    % Default x if empty
    if isempty(x)
        x = 1:length(y);
    end

    % Get axes & color order
    ax = gca;
    if isempty(c)
        co = get(ax,'ColorOrder');
        nLines = numel(findall(ax,'type','line'));
        c = co(mod(nLines, size(co,1))+1, :);
    end

    % --- plot patch and line ---
    h.fill = fill(ax,[x fliplr(x)], [yLower fliplr(yUpper)], c, ...
                  'FaceAlpha', fa, 'EdgeColor','none', ...
                  'DisplayName', description);
    hold on;
    h.line = plot(ax, x, y, 'LineStyle', ls, 'Color', c, 'LineWidth', lw, ...
                  'HandleVisibility','off'); % hide line in legend

    % --- set x-axis scaling ---
    set(ax, 'XScale', xScale);
end

