function [ new_eigenvectors, new_eigenvalues] = IPCA(eigenvectors, eigenvalues, x_t)

    if norm(x_t) < 1e-10
        new_eigenvectors = eigenvectors;
        new_eigenvalues = eigenvalues;
        return;
    end
    % sort
    [~,P] = sort(diag(eigenvalues), 'descend'); eigenvectors = eigenvectors(:,P); eigenvalues = eigenvalues(P,P);

    d = size(eigenvectors,2);
    x_hat = eigenvectors' * x_t;
    x_orth = x_t - eigenvectors * eigenvectors' * x_t;
    
    Q = [eigenvalues + x_hat * x_hat' , norm(x_orth) * x_hat ; norm(x_orth) * x_hat' ,  norm(x_orth)^2];
    
    [U_, S_] = eig(Q);
    [~,idx] = sort(diag(S_), 'descend');
    S_ = S_(idx,idx); 
    U_ = U_(:, idx);
    
    new_eigenvectors = [eigenvectors, x_orth/norm(x_orth)] * U_;
    new_eigenvectors = new_eigenvectors(:,1:d);
    new_eigenvalues = S_(1:d,1:d);
    
    %new_eigenvectors = normc(new_eigenvectors);
    new_eigenvectors = new_eigenvectors ./ vecnorm(new_eigenvectors);
end

