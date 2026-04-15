% make block-lower toeplitz matrix
function BlkTrilToep = make_blk_tril_toeplitz(A,B,C,D,s)
Markovs = get_Markovs(A,B,C,s-1);
BlkTrilToep = Markovs2BlkToeplitz(Markovs,D,depth=s);
end

% make block-lower toeplitz matrix from Markov parameters
function BlkTrilToep = Markovs2BlkToeplitz(Markovs,D,opts)
arguments
    Markovs cell
    D       double
    opts.depth (1,1) double
    opts.toep_struct = [];
end
Markovs2 = [{zeros(size(D))}, {D}, Markovs(:).'];
if isfield(opts,'depth')
    s = opts.depth;
    toep_struct = toeplitz(2:s+1,[2,ones(1,s-1)]);
elseif ~isempty(opts.toep_struct)
    toep_struct = opts.toep_struct;
else
    error('Provide either number of block-rows or toeplitz structure');
end
BlkTrilToep = cell2mat(Markovs2(toep_struct));
end

% get Markov parameters: {CB, CAB, C(A^2)B, ..., C(A^(s-1))B}
function Markovs = get_Markovs(A,B,C,s)
Markovs = cell(1,s);
Markovs{1} = B;
for k = 2:s
     Markovs{k} = A  *  Markovs{k-1};
end
Markovs = cellfun(@(x) C*x, Markovs, UniformOutput=false);
end