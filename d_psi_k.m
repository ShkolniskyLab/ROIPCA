function [res] = d_psi_k(x, k, zkend, dkend, mu, one_min_sum_zk2)
    res = sum((zkend.^2)./(dkend - x).^2) + (one_min_sum_zk2)/(mu - x)^2;
end

