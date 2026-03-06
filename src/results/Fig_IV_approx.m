function fig = Fig_IV_approx(data_type,fig_dir)
fig_file = fullfile(fig_dir,'processed_data.mat');
load(fig_file);
load(fullfile(fig_dir,sprintf('d%s_settings.mat',data_type)));
[p,f,nu,ny,N] = deal(opts.p,opts.f,opts.nu,opts.ny,opts.N);

FS_Tick   = 9;
FS_Label  = 10;
Fs_Legend = 8;

%% Inline make_fig_m1 function

% Extract and set up optional arguments
m1_data = m1;
data_type_plot = 'Re';
X_all = Re_all;
YScale = 'log';
LegCols = [3 2];
LegsLoc = ["northwest","west"];
fontSize = FS_Label;
FS_Legend_val = Fs_Legend;
FigPos = [50 50 252 400];
Units = 'points';
CrameriColors = 'batlow';
fillAlpha = 0.25;
XScale = 'log';
LineWidth = 2;

% Determine x-axis
switch data_type_plot
    case 'N'
        % N_all = X_all;
        xlab = '$N$';
    case {'Re','p'}
        % N_all = repmat(N,1,numel(X_all));
        if strcmp(data_type_plot,'Re')
            xlab = '$\Sigma_e$';
        else
            xlab = '$p$';
        end        
end

% Set colours and make figure
cCram = crameri(CrameriColors, 9);

set_xlabel   = @() xlabel(xlab,'interpreter','latex','FontSize',fontSize);
set_ylabel = @(UorY) ylabel(['$\frac{\big\|\Delta_j \widetilde{',UorY,'}_{\mathrm{f}}\big\|_\mathrm{F}}{\sqrt{N}}$'],...
    'Interpreter','latex','Rotation',0,'FontSize',fontSize*1.25);

fig = figure('Units',Units,'Position',FigPos);
tl2 = tiledlayout(2,1,"TileSpacing",'tight','Padding','tight');
fig.Units = 'pixels';
ax(1) = nexttile();

% --------------------- plotting for Uf_ivs -------------------------------
% Uf_ivs = fieldnames(m1_data.Uf);
% Uf_ivs = setdiff(Uf_ivs,{'iv2a'}); % exclude case that is logically zero
Uf_ivs = {'iv1','iv2c','iv3c','iv4a','iv5a','iv6c'}; % visible 'groups': 3c, 1, 2c, {4a, 4c, 5a, 5c}, 6c
Uf_lb_fs = replace(Uf_ivs,'iv','l');
Uf_ub_fs = replace(Uf_ivs,'iv','u');
Uf_av_fs = replace(Uf_ivs,'iv','m');
num_Uf_ivs = numel(Uf_ivs);

% sqrt_Nall = diag(sqrt(N_all));
u_iv = struct;
for kIV = 1:num_Uf_ivs
    lb_name = Uf_lb_fs{kIV};
    ub_name = Uf_ub_fs{kIV};
    av_name = Uf_av_fs{kIV};
    iv_name = Uf_ivs{kIV};
    
    % u_iv.(av_name) = sqrt_Nall\m1_data.Uf.(iv_name).median;
    % u_iv.(lb_name) = sqrt_Nall\squeeze(m1_data.Uf.(iv_name).pctiles(:,6));
    % u_iv.(ub_name) = sqrt_Nall\squeeze(m1_data.Uf.(iv_name).pctiles(:,16));
    u_iv.(av_name) = m1_data.Uf.(iv_name).median;
    u_iv.(lb_name) = squeeze(m1_data.Uf.(iv_name).pctiles(:,6));
    u_iv.(ub_name) = squeeze(m1_data.Uf.(iv_name).pctiles(:,16));

    % determine color
    switch iv_name
        case 'iv6c', label = '$\hat{\widetilde{U}}_{\mathrm{f,6c}}$'; col = [0 0 0];       LineStyle = '-';
        case 'iv2c', label = '$\hat{\widetilde{U}}_{\mathrm{f,2c}}$'; col = cCram(2,:);    LineStyle = '-.';
        case 'iv3c', label = '$\hat{\widetilde{U}}_{\mathrm{f,3c}}$'; col = cCram(3,:);    LineStyle = '-';
        case 'iv4a', label = '$\hat{\widetilde{U}}_{\mathrm{f,4a}}$'; col = cCram(4,:);    LineStyle = '-*';
        case 'iv4c', label = '$\hat{\widetilde{U}}_{\mathrm{f,4c}}$'; col = cCram(5,:);    LineStyle = '-^';
        case 'iv5a', label = '$\hat{\widetilde{U}}_{\mathrm{f,5a}}$'; col = cCram(6,:);    LineStyle = '-.*';
        case 'iv5c', label = '$\hat{\widetilde{U}}_{\mathrm{f,5c}}$'; col = cCram(7,:);    LineStyle = '-.^';
        case 'iv1',  label = '$U_{\mathrm{f}}$';                  col = cCram(9,:);    LineStyle = '-';
        otherwise,   label = iv_name;                             col = [0 0 0];       LineStyle = '-';
    end

    plotLineWithFill(u_iv.(av_name),u_iv.(lb_name),u_iv.(ub_name), label, ...
        Color=col, FaceAlpha = fillAlpha, LineStyle=LineStyle, LineWidth=LineWidth,...
         x=X_all, XScale=XScale,YScale=YScale);
    hold on;

    u_entries(kIV) = struct('Color',col, 'Alpha',fillAlpha, ...
                    'LineStyle',LineStyle, 'LineWidth',LineWidth, ...
                    'Text',label);
end
set_ylabel('U');
grid on;
axLeg(1) = customLegend(u_entries,ax(1),Location=LegsLoc(1),cols=LegCols(1),FontSize=FS_Legend_val);

% --------------------- plotting for Yf_ivs -------------------------------
ax(2) = nexttile(tl2,2);

Yf_ivs = fieldnames(m1_data.Yf);
Yf_ivs = setdiff(Yf_ivs,{'iv2b'}); % exclude case that is logically zero
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
    
    % y_iv.(av_name) = sqrt_Nall\m1_data.Yf.(iv_name).median;
    % y_iv.(lb_name) = sqrt_Nall\squeeze(m1_data.Yf.(iv_name).pctiles(:,6));
    % y_iv.(ub_name) = sqrt_Nall\squeeze(m1_data.Yf.(iv_name).pctiles(:,16));
    y_iv.(av_name) = m1_data.Yf.(iv_name).median;
    y_iv.(lb_name) = squeeze(m1_data.Yf.(iv_name).pctiles(:,6));
    y_iv.(ub_name) = squeeze(m1_data.Yf.(iv_name).pctiles(:,16));

    % determine color
    switch iv_name
        case 'iv4b', label = '$\hat{\widetilde{Y}}_{\mathrm{f,4b}}$'; col = cCram(4,:);  LineStyle = '-';
        case 'iv5b', label = '$\hat{\widetilde{Y}}_{\mathrm{f,5b}}$'; col = cCram(6,:);  LineStyle = '-.d';
        case 'iv3a', label = '$\Xi_f$' ;                          col = cCram(8,:);  LineStyle = '-x';
        case 'iv6a', label = '$W_f$';                             col = [0 0 0];     LineStyle = '-';
        otherwise,   label = iv_name;                             col = [0 0 0];     LineStyle = '-';
    end

    plotLineWithFill(y_iv.(av_name),y_iv.(lb_name),y_iv.(ub_name), label, ...
        Color=col, FaceAlpha = 0.25, LineStyle=LineStyle, LineWidth=2,...
         x=X_all, XScale=XScale,YScale=YScale);
    hold on;

    y_entries(kIV) = struct('Color',col, 'Alpha',0.25, ...
                    'LineStyle',LineStyle, 'LineWidth',2, ...
                    'Text',label);
end
set_ylabel('Y');
set_xlabel();
grid on;
axLeg(2) = customLegend(y_entries,ax(2),Location=LegsLoc(2),cols=LegCols(2),FontSize=FS_Legend_val);

% Post-processing
% ax(1).YLim(2) = 2;
fig.Units = 'points';

for k=1:2
    ax(k).FontSize = FS_Tick;
    ax(k).YLabel.FontSize = FS_Label*1.5;
end
ax(k).XLabel.FontSize = FS_Label*1.5;

end

%% Helper functions
function color_picker(cCram,IV_name)
    num_colors = size(cCram,1);
    switch IV_name(3:end)
        case '1'
            col = cCram(1,:);
        case {'2a','2b'}
            col = cCram(2,:);
        case '2c'
            col = cCram(3,:);
        case '3a'
            col = cCram(4,:);
        case '3c'
            col = cCram(5,:);
        case {'4a','4b'}
            col = cCram(6,:);
        case '4c'
            col = cCram(7,:);
        case '5c'
            col = cCram(8,:);
        case {'6a','6c'}
            col = zeros(1,3);
    end
end