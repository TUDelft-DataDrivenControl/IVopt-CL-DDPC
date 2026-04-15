function BlkToep = blk_toeplitz_mean(BlkToep,n1,n2)
% averages blocks of size (n1 x n2) over their respective block-diagonals

% determine number of block rows and columns
[S1,S2] = size(BlkToep);
s1 = S1/n1;
s2 = S2/n2;
diaglen_max = min(s1,s2); % maximum block-diagonal length

% iterate over diagonals
for kd = -s1+2:s2-2 % not [-s1+1:s2-1] b/c no averaging occurs for block-diags of length 1
    % determine # blocks on diagonal
    if kd <= 0
        si = -kd;
        diaglen = min(s1-si,diaglen_max);
        i_base = cumsum(n1*ones(1,diaglen))+(si-1)*n1;
        j_base = cumsum(n2*ones(1,diaglen))-n2;
    else
        sj = kd;
        diaglen = min(s2-sj,diaglen_max);
        i_base = cumsum(n1*ones(1,diaglen))-n1;
        j_base = cumsum(n2*ones(1,diaglen))+(sj-1)*n2;
    end
    
    % iterate over entries within block-diagonal (by position within block)
    for ki = 1:n1
        idx = i_base + ki;     % row-indices
        for kj = 1:n2
            jdx = j_base + kj; % column-indices
            
            ind = sub2ind([S1,S2],idx,jdx);                   % get indices of elements of block-diagonal
            A = BlkToep(ind);                          % get entries of block-diagonal
            BlkToep(ind) = repmat(mean(A),diaglen,1);  % reassign original matrix averaged values
        end
    end
end
    
end


