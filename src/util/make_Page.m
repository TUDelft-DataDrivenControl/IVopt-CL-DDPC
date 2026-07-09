function Page = make_Page(data, s1, s2)
% Efficiently constructs Page matrices from input data of the form
% | u_1    u_{1+s2}  ... u_{1+s2*(N-1)}  |
% | ...      ...     ...      ...        |
% | u_{s1} u_{s1+s2} ... u_{s1+s2*(N-1)}a =  |
% 
% where N is the number of columns
%
% Inputs:
%   data: (ndata x Nbar) time series data matrix
%   s1: number of block rows
%   s2: number of sampes for column shift
% Outputs:
%   Page: (ndata*s1 x (floor((Nbar-s1)/s2)+1) full Page matrix

[ndata, Nbar] = size(data);
N = floor((Nbar-s1)/s2)+1; % determine number of columns
Nbar2 = s1+s2*(N-1);       % determine max usable data samples
if Nbar2 ~= Nbar
    error(['Number of provided samples (%d) does not correspond with the number' ...
        'of samples employed by a Page matrix with the specified dimensions (%d)'],Nbar,Nbar2)
end

data = data(:);
Page = zeros(s1*ndata,N);
idxs = [1 s1*ndata];
for k = 1:N
    Page(:,k) = data(idxs(1):idxs(2),1);
    idxs = idxs + s2*ndata;
end
end