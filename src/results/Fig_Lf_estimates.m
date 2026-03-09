function fig = Fig_Lf_estimates(Cases,mLf,iX,opts)
FS_Tick  = 8;
FS_Label = 9;

nCases = numel(Cases);
iXstr = sprintf('iX%d',iX);

Lf_act = mLf.actLf.(iXstr);
for k = 1:nCases
    switch Cases{k}
        case 'actLf'
            Lf_plt = Lf_act;
        otherwise
            Lf_plt = mLf.(Cases{k}).(iXstr).mean - Lf_act;
    end
    maxabsval = max(abs(Lf_plt),[],'all');
end

fig = figure();
fig.OuterPosition(4) = ceil(108.8333*nCases);
tl = tiledlayout(nCases,1,"TileSpacing",'tight','Padding','tight');
for k = 1:nCases
    switch Cases{k}
        case 'actLf'
            Lf_plt = Lf_act;
            yLabel = 'actual';
        otherwise
            Lf_plt = mLf.(Cases{k}).(iXstr).mean - Lf_act;
            if strcmp(Cases{k},'TrPred')
                yLabel = 'TP';
            elseif strcmp(Cases{k},'CLSPC')
                yLabel = 'CL-SPC';
            else
                yLabel = ['IV',Cases{k}(3:end)];
            end
    end
    ax = nexttile;
    imagesc(Lf_plt);
    clim([-maxabsval maxabsval]);
    cmap = crameri('vik','pivot',0);
    colormap(cmap);
    xline(opts.p*opts.nu+0.5); xline(opts.p*(opts.nu+opts.ny)+0.5)
    ylabel(yLabel,'Interpreter','latex','FontSize',FS_Label);
    ax.FontSize = FS_Tick;
    if k~=nCases
        ax.XTickLabel = {};
    end
end
cb = colorbar;
cb.FontSize = FS_Tick;
cb.Layout.Tile = 'east';

fig.Units = 'points';
fig.OuterPosition(3:4) = [252 252/fig.OuterPosition(3)*fig.OuterPosition(4)];
drawnow;

cb.Location = 'manual';
cb.Units = 'points';
tl.Units = 'points';
cbPos = cb.Position;
cb.Position([1 3]) = [cbPos(1)+cbPos(3)/2, cbPos(3)/2]; drawnow;
cbSpace = 5;
tl.OuterPosition(3) = cb.Position(1) - cbSpace - tl.OuterPosition(1);

end

