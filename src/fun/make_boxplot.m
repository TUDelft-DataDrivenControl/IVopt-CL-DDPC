function make_boxplot(data_dir,opts)
arguments (Input)
    data_dir char
    opts.FontSize (1,1) double = 12;
end
% creates nice box plots based on seed_%d.mat files in specified directory

% ------------------------- find .mat files -------------------------------
cwd = pwd;    % save location
cd(data_dir); % move to new location
matfiles = dir('seed_*.mat'); % find seed_%d.mat files
cd(cwd);      % move back to old location

numfiles = length(matfiles);

% names of the controllers
CzNames = {'$\mathcal{Z}_\mathrm{ol}$';...                  1) open-loop IV
           '$\mathcal{Z}^*$';...                            2) optimal IV
           '$\mathcal{Z}_\mathrm{lcf}$';...                 3) IV by LCF
           '$\widehat{\mathcal{Z}}^*_\mathrm{nc}$';...      4) approx. opt. IV w/o controller info
           '$\widehat{\mathcal{Z}}^*_\mathrm{c,lcf}$';...   5) approx. opt. IV w/ controller info - init. w/ 3)
           'CL-SPC';...                                     6) CL-SPC
           '$L_f^*$';...                                    7) using actual L_f
           '$\mathcal{Z}_\mathrm{b}$';...                   8) basic IV
           '$\widehat{\mathcal{Z}}^*_\mathrm{c,b}$'}; %     9) approx. opt. IV w/ controller info - init. w/ 8)

% ------------------------- extracting the data ---------------------------
cost_tot = nan(9,numfiles);
for kf = 1:numfiles
    fp = fullfile(data_dir,matfiles(kf).name); % file path
    cost_tot(:,kf) = load(fp).cost_tot;     % load data
    % cost_tot(:,kf) = load(fp).FroIDerror(:,2);
end

% figure(5);
boxplot(cost_tot(:),repmat(CzNames(:),numfiles,1));
ax = gca;
ax.FontSize = opts.FontSize;
ax.TickLabelInterpreter = 'latex';
% ylim([0 1000]);
ylabel('$\bar{\mathcal{J}}$','interpreter','latex');
grid on;
end