% produce colour-scaled imagesc using vik
function imagesc_vik(CorMat,varargin)
narginchk(1,3)
if nargin >= 2
    vik = varargin{1};
else
    load('vik.mat');
    axh = gca;
end
if nargin >= 3
    axh = varargin{2};
else
    axh = gca;
end
imagesc(axh,CorMat);
c_absmax = max(abs(CorMat),[],'all');
colormap(vik);
colorbar;
caxis([-c_absmax c_absmax]);
end