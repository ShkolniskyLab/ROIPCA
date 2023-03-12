function [new_eigenvectors] = hebbian(eigenvectors, v, itr, eta)

    [n,m] = size(eigenvectors);
    new_eigenvectors = eigenvectors;
    learning_rate = eta/itr;
    
    for i = 1:m
        term = (v' * eigenvectors(:, i)) * v;
        for j = 1:i
                term = term - (v' * eigenvectors(:, i)) * (v' * eigenvectors(:, j)) *  eigenvectors(:, j);
        end
        term = learning_rate * term;
        new_eigenvectors(:, i) = new_eigenvectors(:, i) + term;
    end

    d = size(eigenvectors,1);
    if mod(itr, d) == 0
        new_eigenvectors = orth(new_eigenvectors);
    else
        new_eigenvectors = normc(new_eigenvectors);
    end
    
end

