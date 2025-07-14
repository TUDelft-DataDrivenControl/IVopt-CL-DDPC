function calc_var_mean_est(opts)
arguments (Input)
    opts.FontSize (1,1) double = 12;
end
% creates nice box plots based on seed_%d.mat files in specified directory

% ------------------------- find .mat files -------------------------------
matfiles = dir('seed_*.mat'); % find seed_%d.mat files

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
for kf = 1:numfiles
    fp = fullfile(pwd,matfiles(kf).name);
    Lfs = load(fp).Lfs;
    if kf ==1
        nLf = numel(Lfs{1});
        [s1,s2] = size(Lfs{1});
        LFs = nan(9,nLf,numfiles);
        p = load(fp).p;
        plant = load(fp).plant;
        ny = size(plant.C,1);
        nu = size(plant.B,2)-ny;
    end
    LFs(:,:,kf) = cell2mat(cellfun(@(x) x(:).',Lfs,'UniformOutput',false));
end

% calculate variance
VarLfs = var(LFs,[],3);
var_cmax = max(VarLfs,[],"all");
VarLfs = reshape(VarLfs,9,s1,s2);

% calculate mean error
MeanErrorLfs = mean(LFs-LFs(7,:,:),3);
mean_cmax = max(abs(MeanErrorLfs),[],'all');
MeanErrorLfs = reshape(MeanErrorLfs,9,s1,s2);

% Create tiled layout
tiledlayout(9, 2, 'TileSpacing', 'compact', 'Padding', 'loose');

axVar = gobjects(9,1);
axErr = gobjects(9,1);

text_opts = {'Units', 'normalized', 'Color', 'k',...
    'FontSize', opts.FontSize-3, 'HorizontalAlignment', 'right', ...
    'VerticalAlignment', 'top'};

for kf = 1:9
    % LEFT COLUMN (VarLfs)
    VarLf = squeeze(VarLfs(kf,:,:));
    normVarLf1 = norm(VarLf(:,1:p*nu),"fro");
    normVarLf2 = norm(VarLf(:,p*nu+(1:p*ny)),"fro");
    normVarLf3 = norm(VarLf(:,p*(nu+ny)+1:end),"fro");
    axVar(kf) = nexttile;
    imagesc(VarLf);
    clim([0 var_cmax]); 
    cmap1 = crameri('-devon'); colormap(axVar(kf), cmap1);
    xline(p*nu+0.5); xline(p*(nu+ny)+0.5);
    ylabel(CzNames{kf}, 'Interpreter','latex','FontSize', opts.FontSize);
    text(axVar(kf),0.33, 0.91,sprintf('%.1e',normVarLf1),text_opts{:});
    text(axVar(kf),0.665,0.91,sprintf('%.1e',normVarLf2),text_opts{:});
    text(axVar(kf),0.995,0.91,sprintf('%.1e',normVarLf3),text_opts{:});

    % RIGHT COLUMN (MeanErrorLfs)
    MeanErrorLf = squeeze(MeanErrorLfs(kf,:,:));
    normMeanErLf1 = norm(MeanErrorLf(:,1:p*nu),"fro");
    normMeanErLf2 = norm(MeanErrorLf(:,p*nu+(1:p*ny)),"fro");
    normMeanErLf3 = norm(MeanErrorLf(:,p*(nu+ny)+1:end),"fro");
    axErr(kf) = nexttile;
    imagesc(MeanErrorLf);
    clim([0 mean_cmax]); 
    cmap2 = crameri('-devon'); colormap(axErr(kf), cmap2);
    xline(p*nu+0.5); xline(p*(nu+ny)+0.5);
    text(axErr(kf),0.33, 0.91,sprintf('%.1e',normMeanErLf1),text_opts{:});
    text(axErr(kf),0.665,0.91,sprintf('%.1e',normMeanErLf2),text_opts{:});
    text(axErr(kf),0.995,0.91,sprintf('%.1e',normMeanErLf3),text_opts{:});
end

% Create shared colorbar for VarLfs (left column)
colorbar(axVar(end), 'Position', [0.485 0.11 0.01 0.77]);
title(axVar(1),'Variance','Interpreter','latex','FontSize',opts.FontSize);

% Create shared colorbar for MeanErrorLfs (right column)
colorbar(axErr(end), 'Position', [0.935 0.11 0.01 0.77]);
title(axErr(1),'Mean error','Interpreter','latex','FontSize',opts.FontSize);


end