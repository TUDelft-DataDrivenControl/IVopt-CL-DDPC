function idxs = get_subind_diags(nblkrows,nblkcols,opts)
arguments
    nblkrows (1,1) double {mustBeInteger} % number of block-rows
    nblkcols (1,1) double {mustBeInteger} % number of block-columns
    opts.anti (1,1) logical = false;
    opts.nr (1,1) double {mustBeInteger,mustBeGreaterThanOrEqual(opts.nr,1)} = 1 % number of rows per block-row
    opts.nc (1,1) double {mustBeInteger,mustBeGreaterThanOrEqual(opts.nc,1)} = 1 % number of columns per block-column
end
nr = opts.nr;
nc = opts.nc;

ndiags = nblkrows+nblkcols-1;
selmat = toeplitz(nblkrows:-1:1,nblkrows:ndiags);
if opts.anti % applied to get indices for anti-diagonals
    selmat = flipud(selmat);
end

nrows = nblkrows*nr;
ncols = nblkcols*nc;
ROWs = repmat((1:nrows).',1,ncols);
COLs = repmat( 1:ncols, nrows, 1);

% account for possible multi-dimensionality
if nr > 1 || nc > 1
    selmat = kron(selmat,ones(nr,nc));
end

idxs = cell(1,ndiags);
for ii = 1:ndiags
    row_idxs = ROWs(selmat==ii);
    col_idxs = COLs(selmat==ii);
    idxs{ii} = [row_idxs(:) col_idxs(:)];
end

end