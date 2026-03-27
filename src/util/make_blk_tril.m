function BlkTrilMat = make_blk_tril(Mat,blk_size,varargin)
% turns provided matrix into a block lower triangular form, where
%   Mat         = matrix to block-lower-triangularize
%   blk_size    = size of blocks [s1, s2]
%   varargin{1} = diagonal index to use with tril (optional)

arguments
    Mat double {ismatrix}
    blk_size (1,2) double
end
arguments (Repeating)
    varargin (1,1) double {mustBeInteger}
end
narginchk(2,3)
if nargin == 3
    ndiag = varargin{1}; % diagonal number
else
    ndiag = 0;
end
[m1,m2] = size(Mat);

% block dimensions
b1 = blk_size(1); 
b2 = blk_size(2);

% number of blocks per dimension
s1 = m1/b1;
s2 = m2/b2;
try
    validateattributes([s1 s2],{'double'},{'integer'});
catch
    error('Block dimensions incompatible with provided matrix');
end

BlkTrilMat = Mat.*kron( tril(ones(s1,s2),ndiag), ones(b1,b2) );

end