function make_fig_m1(m1,N_all,opts)
arguments
    m1 struct
    N_all double {mustBeVector}
    opts.Re_all double = [];
    opts.fontSize (1,1) double = 15;
    opts.CrameriColors char = 'roma';  % possible colors: see the 'crameri' command
    opts.FigPos (1,4) double {mustBeReal,mustBeFinite} = [2600 500 1000 600]; % figure position
    opts.fillAlpha (1,1) double {mustBeGreaterThanOrEqual(opts.fillAlpha,0),...
                                 mustBeLessThanOrEqual(opts.fillAlpha,1)} = 0.25
    opts.XScale char {mustBeMember(opts.XScale,["log","linear"])} = 'log';
    opts.YScale char {mustBeMember(opts.YScale,["log","linear"])} = 'linear';
    opts.LineWidth (1,1) double {mustBeFinite,mustBeReal,mustBeGreaterThan(opts.LineWidth,0)} = 2;
end

fontSize  = opts.fontSize;      % font size of x & y labels (scaled later for y due to orientation)
cColors   = opts.CrameriColors;
Position  = opts.FigPos;
fillAlpha = opts.fillAlpha;
XScale    = opts.XScale;
YScale    = opts.YScale;
LineWidth = opts.LineWidth;
Re_all    = opts.Re_all;

% Determine x-axis: N_all (default) or Re_all (if provided)
if ~isempty(Re_all)
    if ~isvector(Re_all)
        error('Re_all must be a vector if not left empty.');
    end
    if numel(unique(N_all)) ~= 1
        error('If Re_all is provided, N_all must have a single unique value.');
    end
    N_all = repmat(N_all,1,numel(Re_all));
    xvals = Re_all;
    xlab = '$\mathrm{Var}(e_k)$';
else
    xvals = N_all;
    xlab = '$N$';
end

% --------------------- set colours and make figure -----------------------
try
    cCram = crameri(cColors, 9); % colors
catch ME
    msg = ME.message;
    msg = sprintf('Crameri color palette not recognized, resulting in:\n%s',msg);
    ME.message = msg;
   rethrow(ME)
end


set_xlabel_2   = @() xlabel(xlab,'interpreter','latex','FontSize',fontSize);
set_ylabel_2_1 = @() ylabel('$\frac{\|\Delta_j U_{\mathrm{f}}^{\mathrm{iv},2a}\|_\mathrm{F}}{\sqrt{N}}$',...
    'Interpreter','latex','Rotation',0,'FontSize',fontSize*1.25);
set_ylabel_2_2 = @() ylabel('$\frac{\|\Delta_j Y_{\mathrm{f}}^{\mathrm{iv},2b}\|_\mathrm{F}}{\sqrt{N}}$',...
    'Interpreter','latex','Rotation',0,'FontSize',fontSize*1.25);

fig2 = figure();
tl2 = tiledlayout(2,1,"TileSpacing",'compact','Padding','compact');
fig2.Units = 'pixels';
fig2.Position = Position;
ax2_1 = nexttile();

% --------------------- plotting for Uf_ivs -------------------------------
% plotting Uf frobenius norm errors
% select the 25th and 75th percentile -> 6 & 16th of 21
% see size(m1.Uf.iv1.pctiles,2);

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
        case 'iv6c', label = '$j=6c$: ref. + 2SLS';        col = [0 0 0];       LineStyle = '-';
        case 'iv2a', label = '$j=2a$: opt. IV';            col = cCram(1,:);    LineStyle = '-.';
        case 'iv2c', label = '$j=2c$: opt. IV + 2SLS';     col = cCram(2,:);    LineStyle = '-.';
        case 'iv3c', label = '$j=3c$: LCF + 2SLS';         col = cCram(3,:);    LineStyle = '--';
        case 'iv4a', label = '$j=4a$: w/o Cz info';        col = cCram(4,:);    LineStyle = ':';
        case 'iv4c', label = '$j=4c$: w/o Cz info + 2SLS'; col = cCram(5,:);    LineStyle = ':';
        case 'iv5a', label = '$j=5a$: w/ Cz info';         col = cCram(6,:);    LineStyle = '-.';
        case 'iv5c', label = '$j=5c$: w/ Cz info + 2SLS';  col = cCram(7,:);    LineStyle = '-.';
        case 'iv1',  label = '$j=1$: OL-IV';               col = cCram(9,:);    LineStyle = '-.';
        otherwise,   label = iv_name;                      col = [0 0 0];       LineStyle = '-';
    end

    plotLineWithFill(u_iv.(av_name),u_iv.(lb_name),u_iv.(ub_name), label, ...
        Color=col, FaceAlpha = fillAlpha, LineStyle=LineStyle, LineWidth=LineWidth,...
         x=xvals, XScale=XScale,YScale=YScale);
    hold on;

    u_entries(kIV) = struct('Color',col, 'Alpha',fillAlpha, ...
                    'LineStyle',LineStyle, 'LineWidth',LineWidth, ...
                    'Text',label);
end
set_ylabel_2_1();
grid on;
axLeg_2_1 = customLegend(u_entries,ax2_1,Location='east',cols=3);
axLeg_2_1.Position([1,2]) = axLeg_2_1.Position([1,2]) + [35 10];

% --------------------- plotting for Yf_ivs -------------------------------
ax2_2 = nexttile(tl2,2);

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
        case 'iv6a', label = '$j=3a,6a$: ref.';      col = [0 0 0];     LineStyle = '-';
        otherwise,   label = field;              col = [0 0 0];     LineStyle = '-';
    end

    plotLineWithFill(y_iv.(av_name),y_iv.(lb_name),y_iv.(ub_name), label, ...
        Color=col, FaceAlpha = 0.25, LineStyle=LineStyle, LineWidth=2,...
         x=xvals, XScale=XScale,YScale=YScale);
    hold on;

    y_entries(kIV) = struct('Color',col, 'Alpha',0.25, ...
                    'LineStyle',LineStyle, 'LineWidth',2, ...
                    'Text',label);
end
set_ylabel_2_2();
set_xlabel_2();
grid on;
axLeg_2_2 = customLegend(y_entries,ax2_2,cols=3,Location='east');
axLeg_2_2.Position([1,2]) = axLeg_2_2.Position([1,2]) + [35 10];
end
