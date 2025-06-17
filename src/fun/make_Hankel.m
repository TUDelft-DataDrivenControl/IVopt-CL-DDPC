function [Hpf, Hp, Hf] = make_Hankel(data, s1, s2)
% Efficiently constructs past-future Hankel matrices from input data
% Inputs:
%   data: (ndata x Nbar) time series data matrix
%   s1: number of past block rows
%   s2: number of future block rows
% Outputs:
%   Hpf: (ndata*(s1+s2) x (Nbar - s1 - s2 + 1)) full Hankel matrix
%   Hp:  (ndata*s1 x num_cols) past block
%   Hf:  (ndata*s2 x num_cols) future block

[ndata, Nbar] = size(data);
num_cols = Nbar - s1 - s2 + 1;

Hpf = zeros(ndata * (s1 + s2), num_cols);

for i = 1:(s1 + s2)
    Hpf((i-1)*ndata+1:i*ndata, :) = data(:, i:i+num_cols-1);
end

Hp = Hpf(1:ndata*s1, :);
Hf = Hpf(ndata*s1+1:end, :);
end