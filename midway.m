function [root, n_itr] = midway(f_mod, d_f_mod, dd, which_eig, z, lambda, mu, one_min_sum_z2)
                rho_inv = 1/lambda;
                y = (dd(which_eig) + dd(which_eig - 1))/2; % initial guess
                fy = f_mod(y, z, dd, rho_inv, mu, one_min_sum_z2);
                zk = z(1:which_eig-1);
                dk = dd(1:which_eig-1);
                zkend = z(which_eig:end);
                dkend = dd(which_eig:end);
                one_min_sum_zk2 = 1 - sum(z(1:which_eig-1).^2);
                %phi_k = @(x, k, z, d) ( sum((z(1:k-1).^2)./(d(1:k-1) - x))  );
                %d_phi_k = @(x, k, zk, dk) ( sum((zk.^2)./(dk - x).^2)  );
                %psi_k = @(x, k, z, d, mu) ( sum((z(k:end).^2)./(d(k:end) - x)) + (1 - sum(z(1:k-1).^2))/(mu - x) );
                %d_psi_k = @(x, k, zkend, dkend, mu, one_min_sum_zk2) ( sum((zkend.^2)./(dkend - x).^2) + (one_min_sum_zk2)/(mu - x)^2 );
                n_itr = 0;
                
                while abs(fy) > 1e-10 && n_itr < 20
                    n_itr = n_itr + 1;

                    %psi_delta = d(which_eig);
                    %psi_r = psi_k(y, which_eig, z, dd, mu) - (psi_delta - y) * d_psi_k(y, which_eig, z, dd, mu);
                    %psi_s = (psi_delta - y)^2 * d_psi_k(y, which_eig, z, dd, mu);

                    %phi_delta = d(which_eig - 1);
                    %phi_r = phi_k(y, which_eig, z, dd) - (phi_delta - y) * d_phi_k(y, which_eig, z, dd);
                    %phi_s = (phi_delta - y)^2 * d_phi_k(y, which_eig, z, dd);

                    Dk = dd(which_eig) - y;
                    Dk1 = dd(which_eig - 1) - y;

                    fy = f_mod(y, z, dd, rho_inv, mu, one_min_sum_z2);
                    a = (Dk + Dk1) * fy - Dk * Dk1 * d_f_mod(y, z, dd, lambda, mu, one_min_sum_z2);
                    b = Dk * Dk1 * fy;
                    c = fy - Dk * d_psi_k(y, which_eig, zkend, dkend, mu, one_min_sum_zk2) - Dk1 * d_phi_k(y, which_eig, zk, dk);
                    
                    if a > 0
                        eta = (2*b) / (a + (a^2 - 4*b*c)^0.5);
                    else
                        eta = (a - (a^2 - 4*b*c)^0.5)/(2*c);
                    end
                    y = y + real(eta);
                end
                
                root = y;
                
end

