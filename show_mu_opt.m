%{
This script reproduces Figure 1

%}

clear;
n0 = 100;
d = 100;
m = 1;
n_test = 500;
n_etas = 500;
Q = orth(randn(d, m));
D = diag(sort(rand( m,1), 'descend'));
cov_matrix = Q * D * Q' + 1e-5 * eye(d);
        
X0 = mvnrnd(zeros(1,d),cov_matrix, n0);
X0TX0 = X0'*X0;
[Q0, D0] = eigs(X0TX0,m);

Xtest = mvnrnd(zeros(1,d),cov_matrix, n_test);
etas = logspace(-5,1,n_etas);
hold on;
for i = 1:n_test
    x_test = Xtest(i,:)';
    [Q_real, D_real] = eigs(X0TX0 + x_test * x_test', m);
    
    z = Q0' * x_test;
    eta_opt = (1/z^2) * (1 - (D0 / D_real));
    for j = 1:n_etas
        eta = etas(j);
        
        Q_updated = Q0 + eta * ( z * x_test );
        Q_updated = Q_updated / norm(Q_updated);
        err(i,j) = min(norm(Q_real - Q_updated) , norm(Q_real + Q_updated));
        
    end
    plot(log10(etas / eta_opt) , log10(err(i,:)));
end
xlabel('log(\eta / \eta*)');
ylabel('log(error)');

hold off;