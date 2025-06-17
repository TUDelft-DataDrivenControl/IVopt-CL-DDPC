function [means, std_devs] = diag_stats(A)
% get means and std deviations of diagonals of a matrix A
    if isvector(A)
        means = A(:).';
        std_devs = zeros(size(means));
    else
        % Get the size of the matrix
        [m, n] = size(A);
        num_diags = m + n - 1; % Number of diagonals
        
        % Preallocate arrays
        means = zeros(1, num_diags);
        std_devs = zeros(1, num_diags);
    
        % Compute statistics for each diagonal
        for k = 1:num_diags
            % Extract diagonal elements
            diag_vals = diag(A, n - k);
    
            % Compute mean and standard deviation
            means(k) = mean(diag_vals);
            std_devs(k) = std(diag_vals);
        end
        means = fliplr(means);
        std_devs = fliplr(std_devs);
    end
end
