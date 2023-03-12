function [best_eta, err] = find_eta(algo_type, XX, XXn, idx, n0, n_comp)
    d = size(XX,2);
    cv = 5;
    [GT,~] = eigs(XX'*XX, n_comp);
    etas = 5.^[5:10];
    X0 = XXn(randperm(size(XXn,1)),:); X0 = X0(1:n0, :);
    [Q0, D0] = eigs(X0' * X0, n_comp);
    
    if strcmp(algo_type, 'GROUSE') == 1
        %etas = linspace(0, pi, 10);
        for k = 1:cv
            perm = randperm(size(XXn,1));
            Xn = XXn(perm,:);
            
            for i = 1:size(etas,2)
                eta = etas(i);
                Q = Q0;
                for j = 1:size(Xn,1)
                    update = Xn(j,:)';
                    idx_for_omega = find(idx(perm(j), :));
                    %idx_for_omega
                    %idx_for_omega = 1:d;
                    [Q] = GROUSE(update, Q, j, idx_for_omega, eta/j);
                end
                err(k,i) = norm(Q*Q' - GT*GT');
            end
        end
        mean_err = mean(err,1);

        idx = find(mean_err == min(mean_err)); idx = idx(1);
        best_eta = etas(idx);
        fprintf('lr for GROUSE = %f\n', best_eta);
        
    elseif strcmp(algo_type, 'ROIPCA') == 1
        % etas = floor(linspace(1, d, 20));
        etas = (1:6) * n_comp;
        for k = 1:cv
            Xn = XXn( (k-1)*n0 + 1 : k*n0, :);
            the_trace_Xn = trace(Q0 * D0 * Q0');
            for i = 1:size(etas,2)
                eta = etas(i);
                Q = Q0;
                D = D0;
               
                [~, SS, VV] = svds(Xn, eta); 
                Xn_trunc = SS * VV';
                for j = 1:size(Xn_trunc,1)
                    update = Xn_trunc(j,:)';
                    the_trace_Xn = the_trace_Xn + norm(update)^2;
                    [Q, D, ~] = update_eigenspectrum([],norm(update)^2,update/norm(update),Q,D,'112',the_trace_Xn, j);
                end
                err(k,i) = norm(Q*Q' - GT*GT');
            end
        end
        mean_err = mean(err,1);

        idx = find(mean_err == min(mean_err)); idx = idx(1);
        best_eta = etas(idx);
        fprintf('r* for ROIPCA = %d\n', best_eta);     
        
    elseif strcmp(algo_type, 'ORPCA') == 1
        
        lambda1s = 2.^[-10:5];
        lambda2s = 2.^[-10:5];
        
        for k = 1:cv
            Xn = XXn( (k-1)*n0 + 1 : k*n0, :);

            for i = 1:size(lambda1s,2)
            	for j = 1:size(lambda2s,2)
                    lambda1 = lambda1s(i);
                    lambda2 = lambda2s(j);
                    A = zeros(n_comp,n_comp);
                    B = zeros(d,n_comp);
                    Q = Q0;
                    for p = 1:size(Xn,2)
                        z = Xn(p,:)';
                        [Q, A, B] = ORPCA(Q, z, A, B, p, lambda1, lambda2);
                    end
                    err(k,i,j) = norm(Q*Q' - GT*GT');
                end
            end
        end
        mean_err = squeeze(mean(err,1));
        [~,I] = min(mean_err(:));
        [I_row, I_col] = ind2sub(size(mean_err),I);
        best_eta = [lambda1s(I_row), lambda2s(I_col)];
        fprintf('lambda for ORPCA = %f and %f \n', best_eta(1), best_eta(2));
        
        
    elseif strcmp(algo_type, 'Oja') == 1
        etas = 2.^[-15:5];
        n0 = size(XXn,1) / cv;
        for k = 1:cv
            Xn = XXn( (k-1)*n0 + 1 : k*n0, :);
            for i = 1:size(etas,2)
                eta = etas(i);
                Q = Q0;
                for j = 1:size(Xn,1)
                    update = Xn(j,:)';
                    Q = OJA(Q, update, j, eta);
                end
                err(k,i) = norm(Q*Q' - GT*GT');
            end
        end
        mean_err = mean(err,1);
        
        idx = find(mean_err == min(mean_err)); idx = idx(1);
        best_eta = etas(idx);
        fprintf('eta for Oja = %f\n', best_eta);
        
        
    elseif strcmp(algo_type, 'Hebbian') == 1
        %[Q0, ~] = eigs(XXn(1:100,:)' * XXn(1:100,:), m);
        etas = 2.^[-15:5];
        n0 = size(XXn,1) / cv;
        for k = 1:cv
            Xn = XXn( (k-1)*n0 + 1 : k*n0, :);
            for i = 1:size(etas,2)
                eta = etas(i);
                Q = Q0;
                for j = 1:size(Xn,1)
                    update = Xn(j,:)';
                    Q = hebbian(Q, update, j, eta);
                end
                err(k,i) = norm(Q*Q' - GT*GT');
            end
        end
        mean_err = mean(err,1);
        
        idx = find(mean_err == min(mean_err)); idx = idx(1);
        best_eta = etas(idx);
        fprintf('eta for Hebbian = %f\n', best_eta);
       
    
    elseif strcmp(algo_type, 'DBPCA') == 1
        lambda1s = 2.^[-10:5];
        lambda2s = 2.^[-10:5];
        lambda1s = 2.^[-1:-1];
        lambda2s = 2.^[-1:-1];
        
           %dbpca
        
        for k = 1:cv
            Xn = XXn( (k-1)*n0 + 1 : k*n0, :);

            for i = 1:size(lambda1s,2)
            	for j = 1:size(lambda2s,2)
                    lambda1 = lambda1s(i);
                    lambda2 = lambda2s(j);
                    
                    d0 = 10;
                    c = lambda1;
                    db = lambda2;
                    di = d0 ./ (2 *  (1:(size(Xn,1)/20)).^2);
                    sizes = floor((c / db^2) * log(d ./ di));
                    sizes = cumsum(sizes);
                    
                    Q = Q0;
                    for p = 1:size(Xn,2)
                        %z = Xn(p,:)';
                        Q = DBPCA(p, sizes, Q, Xn);
                    end
                    err(k,i,j) = norm(Q*Q' - GT*GT');
                end
            end
        end
        mean_err = squeeze(mean(err,1));
        [~,I] = min(mean_err(:));
        [I_row, I_col] = ind2sub(size(mean_err),I);
        best_eta = [lambda1s(I_row), lambda2s(I_col)];
        fprintf('lambda for DBPCA = %f and %f \n', best_eta(1), best_eta(2));
        
        
                    d0 = 10;
                    c = best_eta(1);
                    db = best_eta(2);
                    di = d0 ./ (2 *  (1:(size(Xn,1)/20)).^2);
                    sizes = floor((c / db^2) * log(d ./ di));
                    best_eta = cumsum(sizes);
                    
    end
end

