function fig = Fig_IV_approx(data_type,fig_dir)
fig_file = fullfile(fig_dir,'processed_data.mat');
load(fig_file);
load(fullfile(fig_dir,sprintf('d%s_settings.mat',data_type)));
[p,f,nu,ny,N] = deal(opts.p,opts.f,opts.nu,opts.ny,opts.N);

FS_Tick   = 9;
FS_Label  = 10;
Fs_Legend = 8;

[fig,ax,axLeg] = make_fig_m1(m1,'Re',Re_all,N=N, YScale='log', LegCols=[2 2],...
    LegLocations=["northwest","west"],fontSize=FS_Label,FS_Legend=Fs_Legend,...
    FigPos=[50 50 252 400],Units='points');
ax(1).YLim(2) = 2;
fig.Units = 'points';

for k=1:2
    ax(k).FontSize = FS_Tick;
    ax(k).YLabel.FontSize = FS_Label*1.5;
end
ax(k).XLabel.FontSize = FS_Label*1.5;
end