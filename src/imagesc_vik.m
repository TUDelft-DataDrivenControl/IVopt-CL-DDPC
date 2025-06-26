% produce colour-scaled imagesc using vik
function imagesc_vik(CorMat,opts)
arguments (Input)
    CorMat double
    opts.axh matlab.graphics.axis.Axes
    opts.cmax (1,1) double
end

if isfield(opts,'axh')
    axh = opts.axh;
else
    axh = gca;
end
if isfield(opts,'cmax')
    c_absmax = opts.cmax;
else
    c_absmax = max(abs(CorMat),[],'all');
end

imagesc(axh,CorMat);
clim([-c_absmax c_absmax]);
cmap = crameri('vik','pivot',0);
colormap(cmap);
colorbar;
end