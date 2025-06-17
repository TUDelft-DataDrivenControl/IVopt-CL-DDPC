fig5 = figure(5); tl5 = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
ax5_1 = nexttile(tl5); nxL = (p+2*f-1)*(nu+ny);
spy([A,B;C,D]); grid on;
yline(p*nu+0.5)
xline(p*nu+0.5)
yline(p*(nu+ny)+0.5)
xline(p*(nu+ny)+0.5)
yline(p*(nu+ny)+(2*f-1)*nu+0.5)
xline(p*(nu+ny)+(2*f-1)*nu+0.5)
yline(nxL+0.5,'LineWidth',2);
xline(nxL+0.5,'LineWidth',2);
yline(nxL+f*nu+0.5);
ax5_2 = nexttile(tl5);
imagesc_vik([A,B;C,D]);
yline(p*nu+0.5)
xline(p*nu+0.5)
yline(p*(nu+ny)+0.5)
xline(p*(nu+ny)+0.5)
yline(p*(nu+ny)+(2*f-1)*nu+0.5)
xline(p*(nu+ny)+(2*f-1)*nu+0.5)
yline(nxL+0.5,'LineWidth',2);
xline(nxL+0.5,'LineWidth',2);
yline(nxL+f*nu+0.5);
fig5.Position = [220 151 904 419];