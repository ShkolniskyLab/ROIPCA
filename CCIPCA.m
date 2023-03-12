function [ new_eigenvectors, new_eigenvalues] = CCIPCA(n_samples, eigenvectors, eigenvalues, x_t)
	l = 1;
    n = n_samples;
    m = size(eigenvectors,2);
    f = (1.0 + l)/(1.0 + n_samples);
    new_eigenvectors = eigenvectors;
    new_eigenvalues = eigenvalues;
    for i = 1:m
		nrm = norm(x_t);
        v = eigenvectors(:, i);
        new_eigenvectors(:, i) = (1 - f) * eigenvalues(i) * v + f * (v' * x_t) * x_t;
        %new_eigenvectors(:, i) =  n/(n+1) * v + (1/(n+1)) + x_t * (x_t' * v/norm(v));
        nrm = norm(new_eigenvectors(:, i));
        new_eigenvalues(i) = nrm;
        x_t = x_t - (x_t' * v) * v;
    end

    %new_eigenvectors = normc(new_eigenvectors);
    new_eigenvectors = new_eigenvectors ./ vecnorm(new_eigenvectors);
end

