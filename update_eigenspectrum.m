function [algo_vecs, algo_vals, time, mu] = update_eigenspectrum(A, lambda, v, org_vecs, org_vals, algo_type, trace_, itr_num)
%	Update eigenspectrum of symmetric matrix A after rank-one
%	pertrubation, i.e., find the eigenspectrum of A_ = A + lambda * v'v;
%	Using the method described in the paper. If the original
%	eigendecompositin is not provided (org_vecs/org_vals) It will be
%	calulated.
% algo_type = a 2-char string of the form (0/1/2)(1/2) indicating the type of
% equation to solve. Left char indicates TSE - 1 or CTSE - 2. Right char
% indicates TEF - 0, CTEF1 - 1, CEEF2 - 2
 
    time = 0;
    % make matrix
    if size(org_vals, 1) == 1 || size(org_vals, 2) == 1
        org_vals = diag(org_vals);
    end

    mitz_type = 0;
    n0 = 1;
    if mitz_type==1
        
        %stam = (n0-1) * diag(shrinker_wrapper_short(diag(org_vals)/(n0-1), 100, (n0-1), 'O_1', 1));
        org_vals = (n0 - 1) * diag(apply_inverse_ev(diag(org_vals/(n0 - 1)), n0 - 1, size(A,1), 1));        %org_vecs = real_Q;
        org_vecs(:, diag(org_vals) == 0) = [];
        org_vals = diag(org_vals); org_vals(org_vals == 0) = []; org_vals = diag(org_vals);

    end
    
    if strcmpi(algo_type, '110') == 1
        se_type = 'se1';
        vf_type = 'ef1';
        mu_type = 'zero';
    elseif strcmpi(algo_type, '111') == 1
        se_type = 'se1';
        vf_type = 'ef1';
        mu_type = 'opt';
    elseif strcmpi(algo_type, '112') == 1
        se_type = 'se1';
        vf_type = 'ef1';
        mu_type = 'mean';
    elseif strcmpi(algo_type, '120') == 1
        se_type = 'se1';
        vf_type = 'ef2';
        mu_type = 'zero';
    elseif strcmpi(algo_type, '121') == 1
        se_type = 'se1';
        vf_type = 'ef2';
        mu_type = 'opt';
    elseif strcmpi(algo_type, '122') == 1
        se_type = 'se1';
        vf_type = 'ef2';
        mu_type = 'mean';
    elseif strcmpi(algo_type, '210') == 1
        se_type = 'se2';
        vf_type = 'ef1';
        mu_type = 'zero';
    elseif strcmpi(algo_type, '211') == 1
        se_type = 'se2';
        vf_type = 'ef1';
        mu_type = 'opt';
    elseif strcmpi(algo_type, '212') == 1
        se_type = 'se2';
        vf_type = 'ef1';
        mu_type = 'mean';
    elseif strcmpi(algo_type, '220') == 1
        se_type = 'se2';
        vf_type = 'ef2';
        mu_type = 'zero';
    elseif strcmpi(algo_type, '221') == 1
        se_type = 'se2';
        vf_type = 'ef2';
        mu_type = 'opt';
    elseif strcmpi(algo_type, '222') == 1
        se_type = 'se2';
        vf_type = 'ef2';
        mu_type = 'mean';
    elseif strcmpi(algo_type, '310') == 1
        se_type = 'se3';
        vf_type = 'ef1';
        mu_type = 'zero';
    elseif strcmpi(algo_type, '312') == 1
        se_type = 'se3';
        vf_type = 'ef1';
        mu_type = 'mean';
    else
        fprintf('unkown algorithm type. Aborting. \n');
        return;
    end
    
    % remove very small entries of v and normalize
%     v(abs(v) < 1e-12) = 0;
%     v = v/norm(v);
    

    
    if size(org_vecs, 2) ~= size(org_vals,1)

        fprintf('number of input eigevalues should be the same as number of input eigenvectors. Aborting.');
        return;
    end
    
%     if sparse == 0 && issparse(A) == 1
%         A = full(A);
%         v = full(v);
%     end
%     
%     if sparse == 1 && issparse(A) == 0
%         A = sparse(A);
%         v = sparse(v);
%     end
    
    if size(org_vecs, 1) == 1 && size(org_vecs, 2) == 1 % org_eigvecs size tell how many pairs to calculate
        fprintf('old eigenvectors/old eigenvalues were not provided. Calculating the first %d eigenpairs... ', org_vecs);
        [org_vecs, org_vals] = eigs(A, org_vecs);
        fprintf('done.\n');
    end

    % descending order
    %[~,P] = sort(diag(org_vals), 'descend'); org_vecs = org_vecs(:,P); org_vals = org_vals(P,P);
    
     if lambda == 0
        algo_vecs = org_vecs;
        algo_vals = diag(org_vals);
        time = 0;
        mu = 0;
        return
    end
    
    % variables and parameters

    n = size(org_vecs, 1); % dimension
    m = size(org_vecs, 2); % number of known eigenpairs
    
    dd = diag(org_vals);
    algo_vals = zeros(m, 1);
    algo_vecs = zeros(n, m);
   
    n_iterations = [];
   
    z = org_vecs'*v;
    
    one_min_z = 1 - real(z'*conj(z));
    
    if (strcmpi(mu_type, 'opt') == 1 || strcmpi(se_type,'se2') == 1)
        s = (v'*A)*(v - org_vecs*(org_vecs'*v));
    else
        s = 0;
    end
    %fprintf('mu = %s', mu);
    if strcmpi(mu_type, 'opt') == 1
        mu = s / one_min_z;
    elseif strcmpi(mu_type, 'mean') == 1    
        mu = (trace_ - sum(dd))/(n-m);
    elseif strcmpi(mu_type, 'zero') == 1    
        mu = 0;
    end
    
    mu = real(mu);
%     mu
%     if mu==0
%         mu = 0.2 * ((trace(A) - sum(dd))/(n-m)) + 0.8 * org_vals(end);
%     end

%     factor1 = ((1 - norm(z)^2));
%     factor2 = ((1 - norm(z)^2)*(trace(A) - sum(dd)))/(n-m);
%     mu = (trace(A) - sum(dd))/(n-m);
%     mu = dd(end)/2;
%     mu = (v'*A)*((eye(n) - org_vecs*org_vecs')*v);
%     mu = 0;
    
    % calculate eigenvalues
    if strcmpi(se_type,'se1') || m == n
        %f =  inline('1 + lambda * sum(sum((u.^2)./(d - x))) + lambda * one_min_z./ (mu - x) ', 'x', 'u', 'd', 'lambda', 'one_min_z', 's', 'mu');
        %df = inline('lambda * sum(sum((u.^2)./((d - x).^2))) + lambda * one_min_z./ ((mu - x)^2) ', 'x', 'u', 'd', 'lambda', 'one_min_z', 's', 'mu');
        f = @(x, u, d, lambda, one_min_z, s, mu) (1 + lambda * sum(sum((real(u.*conj(u)))./(d - x))) + lambda * one_min_z./ (mu - x));
        df = @(x, u, d, lambda, one_min_z, s, mu) (lambda * sum(sum((real(u.*conj(u)))./((d - x).^2))) + lambda * one_min_z./ ((mu - x)^2));
        
        if mitz_type == 2
            if lambda < 0
                f = @(x, u, d, lambda, one_min_z, s, mu) (1 + lambda * sum(sum((u.^2)./(d - x))) - sum(1./(d - x)) + lambda * one_min_z./ (mu - x) + (100 - size(u,1))./(mu - x));
                df = @(x, u, d, lambda, one_min_z, s, mu) (lambda * sum(sum((u.^2)./((d - x).^2))) -  sum(1./(d - x).^2) + lambda * one_min_z./ ((mu - x)^2)+ (100 - size(u,1))./(mu - x)^2);
            end
        end
        
        
        %f = @(x, u, d, lambda, one_min_z, s, mu) (1 + lambda * sum(sum((u.^2)./(d - (1./x)))) + lambda * one_min_z./ (mu - 1/x));
        %df = @(x, u, d, lambda, one_min_z, s, mu) (-lambda * sum(sum((u.^2)./((d.*x - 1).^2))) - lambda * one_min_z./ ((mu.*x - 1)^2));
    elseif strcmpi(se_type,'se2')
        %f =  inline('1 + lambda * sum(sum((u.^2)./(d - x))) + lambda * one_min_z./ (mu - x) - lambda * (s - mu * one_min_z)/((mu - x)^2)', 'x', 'u', 'd', 'lambda', 'one_min_z', 's', 'mu');
        %df = inline('lambda * sum(sum((u.^2)./((d - x).^2))) + lambda * one_min_z./ ((mu - x)^2) - 2 * lambda * (s - mu * one_min_z)/((mu - x)^3)', 'x', 'u', 'd', 'lambda', 'one_min_z', 's', 'mu');
        f = @(x, u, d, lambda, one_min_z, s, mu) (1 + lambda * sum(sum((u.^2)./(d - x))) + lambda * one_min_z./ (mu - x) - lambda * (s - mu * one_min_z)/((mu - x)^2));
        df = @(x, u, d, lambda, one_min_z, s, mu) (lambda * sum(sum((u.^2)./((d - x).^2))) + lambda * one_min_z./ ((mu - x)^2) - 2 * lambda * (s - mu * one_min_z)/((mu - x)^3));
    elseif strcmpi(se_type,'se3')
        f = @(x, u, d, lambda, one_min_z, s, mu) (1 + lambda * sum(sum((u.^2)./(d - x))) + lambda * one_min_z./ (mu - x) +   1^2 * (sum(1./(d-x)) - sum(1./(mu -x)))      );
        df = @(x, u, d, lambda, one_min_z, s, mu) (lambda * sum(sum((u.^2)./((d - x).^2))) + lambda * one_min_z./ ((mu - x)^2)    +   1^2 * (sum(1./(d-x).^2) - sum(1./(mu -x).^2))         );
%          f = @(x, u, d, lambda, one_min_z, s, mu) (1 + lambda * sum(sum((u.^2)./(d - x))) + lambda * one_min_z./ (mu - x) - (1000 * 1^2)/(n0 - x)   );
%          df = @(x, u, d, lambda, one_min_z, s, mu) (lambda * sum(sum((u.^2)./((d - x).^2))) + lambda * one_min_z./ ((mu - x)^2)  - (1000 * 1^2)/(n0 - x)^2  );
    end
    
    bol = 0;
    
    %fprintf('solving %s using newton bisection method... ', se_type); 
    max_itr = 50;
    prox = 1e-10;
    tic;
    one_min_sum_z2 = 1 - sum(z.^2);
            
            if lambda > 0
                
                [algo_vals(1), n_itr] = nr( f, df, dd(1), dd(1) + 1 * lambda, prox, max_itr ,  z(1:m), dd(1:m), lambda, one_min_z, s, mu);

                for which_eig = 2:m
                    %interval = [1.000000001 * dd(which_eig) , 0.*9999999 * dd(which_eig - 1)];
                    %f(interval(1) , z_, d_, lambda_, one_min_z_, s_, mu_)
                    %f(interval(1) , z_, d_, lambda_, one_min_z_, s_, mu_)
                    %algo_vals(which_eig) = fzero(@(x) f(x, z_, d_, lambda_, one_min_z_, s_, mu_),interval);
                    
                    
                    %[algo_vals(which_eig)] = nr( f, df, dd(which_eig), dd(which_eig - 1), prox, max_itr ,  z(1:m), dd(1:m), lambda, one_min_z, s, mu);
                     rho_inv = 1/lambda;
                     f_mod = @(x, z, d, rho_inv, mu, one_min_sum_z2) (rho_inv + sum((z.^2) ./ (d - x)) + (one_min_sum_z2)/(mu - x) );
                     d_f_mod = @(x, z, d, rho_inv, mu, one_min_sum_z2) ( sum((z.^2) ./ (d - x).^2) + (one_min_sum_z2)/(mu - x)^2 );
                     [algo_vals(which_eig), tmp] = midway(f_mod, d_f_mod, dd, which_eig, z, lambda, mu, one_min_sum_z2);
                end 
 
                
                
            else
                fprintf('lambda must be > 0\n');
            end
    
    %algo_vals0 = algo_vals;
    %algo_vals = eigs(A + lambda * (v * v'), size(org_vecs,2));
    %norm(algo_vals - algo_vals0)/norm(algo_vals)

    %mu = 0;
    diag_org_vals = diag(org_vals);
    for which_eig = 1:m
        temp_diag = 1./(diag_org_vals - algo_vals(which_eig));
%         if max(temp_diag) > 100000
%             %temp_diag
%         end
        if (max(temp_diag) == inf || max(temp_diag) == -inf || sum(isnan(temp_diag)) > 0 || max(temp_diag) > 10000)
            algo_vecs = org_vecs;
            algo_vals = org_vals;
            return;
        else
                algo_vecs(:, which_eig) = org_vecs * (temp_diag.* z);
        end
    end 
    
%    if bol == 1; 
%             did_not_converge
%              org_vals
%              algo_vals
%              algo_vecs 
%              
%              return
%     end
    
    if strcmpi(vf_type, 'none') == 0 && m~=n
       
                %fprintf('correcting [%s]... ', vf_type);
                correction_vector = v - org_vecs*(conj(org_vecs')*v);
                 if strcmpi(vf_type,'ef1') == 1
                     
                     %algo_vecs = algo_vecs + correction_vector.*repmat((1./(mu - algo_vals))', n, 1);
                     
                     %repmat((1./(mu - algo_vals))', n, 1) .* correction_vector;
                     
                     for k = 1:m  
                         
                            algo_vecs(:,k) = algo_vecs(:,k) + correction_vector * (1/(mu - algo_vals(k))); 
                     end
                 elseif strcmpi(vf_type,'ef2')
                     b = A*correction_vector; 
                     %b = (algo_vecs * diag(algo_vals) * algo_vecs' + mu * (eye(n) - algo_vecs * algo_vecs')) * correction_vector;   
                     for k = 1:m  
                            
                            algo_vecs(:,k) = algo_vecs(:,k) + correction_vector * (1/(mu - algo_vals(k))) + mu * (1/(mu - algo_vals(k))^2) * correction_vector - (1/(mu - algo_vals(k))^2) * b;
                     end

                 end
    end

    if mod(itr_num, n) == 0
        algo_vecs = orth(algo_vecs);
    else
        algo_vecs = normc(real(algo_vecs));
    end
    
    %algo_vecs = algo_vecs ./ vecnorm(algo_vecs);
    
    
    if rad2deg(subspace(algo_vecs, org_vecs)) > 50
        algo_vecs = org_vecs;
        algo_vals = org_vals;
        %fprintf('Done.\n');
        
    end
    time = toc;
    %fprintf('Done.\n');
    
    
end