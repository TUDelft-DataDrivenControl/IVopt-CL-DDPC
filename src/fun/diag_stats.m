function [means, std_devs] = diag_stats(A,opts)
arguments
    A double {mustBeMatrix}
    opts.nr (1,1) double {mustBeInteger,mustBePositive} = 1
    opts.anti logical = false; % use anti-diagonals?
end
% get means and std deviations of diagonals of a matrix A
nr = opts.nr; % number of rows per block

if opts.anti
    A = flipud(A);
end

if isvector(A)
    A = A(:);
    
    if rem(numel(A),nr)~=0
        error('Block dimensions do not correspond with vector dimensions');
    end

    means = A.';
    std_devs = zeros(size(means));
else
    % Get the size of the matrix
    [m, n] = size(A);
    if rem(m,nr)~=0
        error('Block dimensions do not correspond with number of matrix rows');
    end
    num_diags = m/nr + n - 1; % Number of diagonals
    
    % Preallocate arrays
    means    = zeros(nr, num_diags);
    std_devs = zeros(nr, num_diags);
    
    % loop over block dimensions
    for knr = 1:nr
        Anr = A(knr:nr:end,:);
        % Compute statistics for each diagonal
        for k = 1:num_diags
            % Extract diagonal elements
            diag_vals = diag(Anr, n - k);
    
            % Compute mean and standard deviation
            means(knr,k)    = mean(diag_vals);
            std_devs(knr,k) = std(diag_vals);
        end
    end
    means    = fliplr(means);
    std_devs = fliplr(std_devs);

    if opts.anti
        means    = flipud(means);
        std_devs = flipud(std_devs);
    end

end
end
