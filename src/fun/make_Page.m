function [Hpf, Hp, Hf] = make_Page(data, s1, s2)
% Efficiently constructs past-future Hankel matrices from input data
% Inputs:
%   data: (ndata x Nbar) time series data matrix
%   s1: number of past block rows
%   s2: number of future block rows
% Outputs:
%   Hpf: (ndata*(s1+s2) x num_cols ) full Page matrix
%   Hp:  (ndata*s1 x num_cols) past block
%   Hf:  (ndata*s2 x num_cols) future block

[ndata, Nbar] = size(data);
num_cols = floor(Nbar/(s1+s2));
Nbar = num_cols*(s1+s2);

Hpf = reshape(data(:,1:Nbar),ndata*(s1+s2),num_cols);
Hp = Hpf(1:ndata*s1, :);
Hf = Hpf(ndata*s1+1:end, :);
end