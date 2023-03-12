function [algo_vecs, algo_vals, time, n_itr, lr] = update_eigenspectrum_fast(A, lambda, v, org_vecs, org_vals, algo_type, trace_, itr_num)
%	Update eigenspectrum of symmetric matrix A after rank-one
%	pertrubation, i.e., find the eigenspectrum of A_ = A + lambda * v'v;
%	Using the method described in the paper. If the original
%	eigendecompositin is not provided (org_vecs/org_vals) It will be
%	calulated.
% algo_type = a 2-char string of the form (0/1/2)(1/2) indicating the type of
% equation to solve. Left char indicates TSE - 1 or CTSE - 2. Right char
% indicates TEF - 0, CTEF1 - 1, CEEF2 - 2
    lr = 1;
    time = 0;
    n_itr = 1;
    if size(org_vals, 1) == 1 || size(org_vals, 2) == 1
        org_vals = diag(org_vals);
    end    
    
    if strcmpi(algo_type, '110') == 1
        se_type = 'se1';
        vf_type = 'ef1';
        mu = 'zero';
    elseif strcmpi(algo_type, '111') == 1
        se_type = 'se1';
        vf_type = 'ef1';
        mu = 'opt';
    elseif strcmpi(algo_type, '112') == 1
        se_type = 'se1';
        vf_type = 'ef1';
        mu = 'mean';
    elseif strcmpi(algo_type, '113') == 1
        se_type = 'se1';
        vf_type = 'ef1';
        mu = 'opt_approx';
    elseif strcmpi(algo_type, '120') == 1
        se_type = 'se1';
        vf_type = 'ef2';
        mu = 'zero';
    elseif strcmpi(algo_type, '121') == 1
        se_type = 'se1';
        vf_type = 'ef2';
        mu = 'opt';
    elseif strcmpi(algo_type, '122') == 1
        se_type = 'se1';
        vf_type = 'ef2';
        mu = 'mean';
    elseif strcmpi(algo_type, '210') == 1
        se_type = 'se2';
        vf_type = 'ef1';
        mu = 'zero';
    elseif strcmpi(algo_type, '211') == 1
        se_type = 'se2';
        vf_type = 'ef1';
        mu = 'opt';
    elseif strcmpi(algo_type, '212') == 1
        se_type = 'se2';
        vf_type = 'ef1';
        mu = 'mean';
    elseif strcmpi(algo_type, '220') == 1
        se_type = 'se2';
        vf_type = 'ef2';
        mu = 'zero';
    elseif strcmpi(algo_type, '221') == 1
        se_type = 'se2';
        vf_type = 'ef2';
        mu = 'opt';
    elseif strcmpi(algo_type, '222') == 1
        se_type = 'se2';
        vf_type = 'ef2';
        mu = 'mean';
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
        return
    end
    
    % variables and parameters

    n = size(org_vecs, 1); % dimension
    m = size(org_vecs, 2); % number of known eigenpairs
    
    dd = diag(org_vals);
    algo_vals = zeros(m, 1);
    algo_vecs = zeros(n, m);
    

    z = org_vecs'*v;
    
    one_min_z = 1 - z'*z;
    if (strcmpi(mu, 'opt') == 1 || strcmpi(se_type,'se2') == 1)
        s = (v'*A)*(v - org_vecs*(org_vecs'*v));
    else
        s = 0;
    end
    %fprintf('mu = %s', mu);
    if strcmpi(mu, 'opt') == 1
        mu = s / one_min_z;
        
    elseif strcmpi(mu, 'mean') == 1    
        mu = (trace_ - sum(dd))/(n-m);
    elseif strcmpi(mu, 'zero') == 1    
        mu = 0;
    end

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
        %f = @(x, u, d, lambda, one_min_z, s, mu1, mu2, i) (1 + lambda * u(i).^2 / (d(i) - x) + lambda * sum(sum((u(setdiff(1:end,i)).^2))./(mu1 - x)) + lambda * one_min_z./ (mu2 - x));
        %df = @(x, u, d, lambda, one_min_z, s, mu1, mu2, i) (lambda * u(i).^2 / (d(i) - x)^2 + lambda * sum(sum((u(setdiff(1:end,i)).^2))./((mu1 - x).^2)) + lambda * one_min_z./ ((mu2 - x)^2));
        f = @(x, u, d, lambda, one_min_z, s, mu) (1 + lambda * sum(sum((u.^2)./(d - x))) + lambda * one_min_z./ (mu- x));
        df = @(x, u, d, lambda, one_min_z, s, mu) (lambda * sum(sum((u.^2)./((d - x).^2))) +  lambda * one_min_z./ ((mu - x)^2));
   
    elseif strcmpi(se_type,'se2')
        %f =  inline('1 + lambda * sum(sum((u.^2)./(d - x))) + lambda * one_min_z./ (mu - x) - lambda * (s - mu * one_min_z)/((mu - x)^2)', 'x', 'u', 'd', 'lambda', 'one_min_z', 's', 'mu');
        %df = inline('lambda * sum(sum((u.^2)./((d - x).^2))) + lambda * one_min_z./ ((mu - x)^2) - 2 * lambda * (s - mu * one_min_z)/((mu - x)^3)', 'x', 'u', 'd', 'lambda', 'one_min_z', 's', 'mu');
        f = @(x, u, d, lambda, one_min_z, s, mu1, mu2, i) (1 + lambda * sum(sum((u.^2)./(d - x))) + lambda * one_min_z./ (mu - x) - lambda * (s - mu * one_min_z)/((mu - x)^2));
        df = @(x, u, d, lambda, one_min_z, s, mu1, mu2, i) (lambda * sum(sum((u.^2)./((d - x).^2))) + lambda * one_min_z./ ((mu - x)^2) - 2 * lambda * (s - mu * one_min_z)/((mu - x)^3));
    end
    
 
    bol = 0;
    max_itr = 100;
    prox = 1e-1;
    n_iterations = [];
    
   	tic;
    %fprintf('solving %s using newton bisection method... ', se_type); 
            %mu_known = ((sum(dd) - dd(1)) / (m -1));
            %mu_known = dd(setdiff(1:end,1));
            %mu_known = dd(1);
            [algo_vals(1)] = nr( f, df, dd(1), dd(1) + 1 * lambda, prox, max_itr ,  z(1:m), dd(1:m), lambda, one_min_z, s,  mu);
           
            one_min_sum_z2 = 1 - sum(z.^2);
            n_itr_midway = [];  n_itr_fast = [];
            for which_eig = 2:m
                rho_inv = 1/lambda;
                f_mod = @(x, z, d, rho_inv, mu, one_min_sum_z2) (rho_inv + sum((z.^2) ./ (d - x)) + (one_min_sum_z2)/(mu - x) );
                d_f_mod = @(x, z, d, rho, mu, one_min_sum_z2) ( sum((z.^2) ./ (d - x).^2) + (one_min_sum_z2)/(mu - x)^2 );
                [algo_vals(which_eig), n_itr] = midway(f_mod, d_f_mod, dd, which_eig, z, lambda, mu, one_min_sum_z2);
                %n_itr_midway = [n_itr_midway, tmp];
                %n_iterations = [n_iterations, n_itr];
                %mu_known = mean(mean(dd(setdiff(1:end,which_eig))));
                %mu_known = dd(setdiff(1:end,which_eig));
                %mu_known = dd(which_eig);
                %[algo_vals(which_eig), tmp] = nr( f, df, dd(which_eig), dd(which_eig - 1), prox, max_itr ,  z(1:m), dd(1:m), lambda, one_min_z, s, mu);
                %n_itr_fast = [n_itr_fast, tmp];
            end 
       
    %fprintf('Done.\n');
    %algo_vals = eigs(A + lambda + v*v', m);
    algo_vals = real(algo_vals);
    % calculate eigenvectors
    %fprintf('calculating eigenvectors... ');
    
    
    diag_org_vals = diag(org_vals);
    %QQTv = org_vecs * z;
    correction_vector2 = v - org_vecs*(org_vecs'*v);
    
    %%
    flag = 0;
    for which_eig = 1:m
        lr(which_eig) = (1./z(which_eig)^2) * (1 - org_vals(which_eig, which_eig) /  algo_vals(which_eig));
        %lr(which_eig) = (1 - org_vals(which_eig, which_eig) /  algo_vals(which_eig));
    end
    
    flag = 0;
    %%
    
    for which_eig = 1:m
        temp_diag = 1./(diag_org_vals - algo_vals(which_eig));
        
        if (max(temp_diag) == inf || max(temp_diag) == -inf || sum(isnan(temp_diag)) > 0)
            algo_vecs = org_vecs;
            algo_vals = org_vals;
            return;
        else

                
                  beta = z(which_eig) *  (1/(org_vals(which_eig, which_eig) - algo_vals(which_eig)));   
                  algo_vecs(:,which_eig) = org_vecs(:,which_eig); 
                  algo_vecs(:,which_eig) = algo_vecs(:,which_eig) + (1/beta) * correction_vector2 * (1/(mu - algo_vals(which_eig))); 
    
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
                %correction_vector2 = v - org_vecs*(org_vecs'*v);
                 if strcmpi(vf_type,'ef1') == 1
                     %algo_vecs = algo_vecs + correction_vector.*repmat((1./(mu - algo_vals))', n, 1);
                     
                     %repmat((1./(mu - algo_vals))', n, 1) .* correction_vector;
                     
                     for k = 1:m  
                            %correction_vector2 = v - org_vecs*(org_vecs'*v); % TMP
                            %mu1 = (sum(sum(dd)) - dd(k))/(m - 1);
                            %mu = (trace(A) - sum(dd) + dd(k))/(n-1);
                            %tau = 1./(diag(org_vals) - algo_vals(k)); tau(k) = 0; tau = sum(tau)/(m-1);
                             tau = 0;
                            % beta = z(k) * (1 / (org_vals(k, k) - algo_vals(k)) - tau);
                            %algo_vecs(:,k) = algo_vecs(:,k) + correction_vector2 * (1/(mu - algo_vals(k))); 
                     
%                             v = v - (v' *  org_vecs(:,k)) * org_vecs(:,k) ;
%                             z = org_vecs'*v;
                     end
                 elseif strcmpi(vf_type,'ef2')
                     
                     %b = A*correction_vector2; 
                     %b = (algo_vecs * diag(algo_vals) * algo_vecs' + mu * (eye(n) - algo_vecs * algo_vecs')) * correction_vector;   
                     for k = 1:m  
                            
                            algo_vecs(:,k) = algo_vecs(:,k) + correction_vector2 * (1/(mu - algo_vals(k))) + mu * (1/(mu - algo_vals(k))^2) * correction_vector2 - (1/(mu - algo_vals(k))^2) * b;
                     end
                 end
    end
    
%     if max(lr) > 100
%         algo_vecs = org_vecs;
%         algo_vals = diag(org_vals);
%     end
    
    if mod(itr_num, n) == 0
        algo_vecs = orth(algo_vecs);
    else
        algo_vecs = normc(real(algo_vecs));
    end
    
%     if rad2deg(subspace(algo_vecs, org_vecs)) > 50
%         algo_vecs = org_vecs;
%         algo_vals = org_vals;
%         
%         %fprintf('hey\n');
%     end
    time = toc;
    %fprintf('Done.\n');
    
    
end