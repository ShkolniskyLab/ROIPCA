function [ new_eigenvectors] = OJA(eigenvectors, x_t, i, eta)

    %x_t = x_t';

    new_eigenvectors = eigenvectors + (eta/i) * x_t * (x_t' * eigenvectors);

    %[new_eigenvectors, ~] = qr(new_eigenvectors);
    
    
    d = size(eigenvectors,1);
    if mod(i, d) == 0
        new_eigenvectors = orth(new_eigenvectors);
    else
        new_eigenvectors = normc(new_eigenvectors);
    end
    
end

