function [res] = d_phi_k(x, k, zk, dk)
    res = sum((zk.^2)./(dk - x).^2);
end

