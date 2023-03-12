%{
This script reproduces Figure 2

%}

n = 1000;
m = 10;
n_for_stats = 10;
ds = [100,200,500,1000,2000,5000, 10000];
summary = zeros(size(ds,2), 4);
n_d = size(ds,2);
for i = 1:n_d
    d = ds(i);
    d
    mu_ = zeros(d,1);
    cov_ = diag( [ 1 + rand(m, 1) ; zeros(d - m, 1)]);
    
    X = mvnrnd(mu_, cov_, n);
    XX = X'*X;
    [Q,D] = eigs(XX, m); D = diag(D);
    
    total_t_roipca1 = 0;
    total_t_roipca2 = 0;
    total_t_roipca_fast = 0;
    total_t_roipca_fast2 = 0;
    total_t_ipca = 0;
    total_t_ccipca = 0;
    for j = 1:n_for_stats
        
        update = mvnrnd(mu_, cov_, 1)';
        update_norm = update/norm(update);
        
        [~,~, t_roipca1] = update_eigenspectrum(XX,norm(update)^2,update_norm,Q,diag(D), '110',[], j);
        %[~,~, t_roipca2] = update_eigenspectrum(XX,norm(update)^2,update_norm,Q,diag(D), '220',[], j);
        [~,~, t_roipca_fast] = update_eigenspectrum_fast(XX,norm(update)^2,update_norm,Q,diag(D), '110',[], j);
        %[~,~, t_roipca_fast2] = update_eigenspectrum_fast(XX,norm(update)^2,update_norm,Q,diag(D), '220',[], j);
        tic; IPCA(Q, diag(D), update); t_ipca = toc;
        tic; CCIPCA(1 + n - 1, Q, diag(D), update); t_ccipca = toc;
        total_t_roipca1 = total_t_roipca1 + t_roipca1;
%        total_t_roipca2 = total_t_roipca2 + t_roipca2;
        total_t_roipca_fast = total_t_roipca_fast + t_roipca_fast;
%        total_t_roipca_fast2 = total_t_roipca_fast2 + t_roipca_fast2;
        total_t_ipca = total_t_ipca + t_ipca;
        total_t_ccipca = total_t_ccipca + t_ccipca;
    end
    
    summary(i,1) = d;
    summary(i,2) = total_t_roipca1 / n_for_stats;
    summary(i,3) = total_t_roipca2 / n_for_stats;
    summary(i,4) = total_t_roipca_fast / n_for_stats;
    summary(i,5) = total_t_roipca_fast2 / n_for_stats;
    summary(i,6) = total_t_ipca / n_for_stats;
    summary(i,7) = total_t_ccipca / n_for_stats;
end

        hold on;
        plot(log(summary(:,1)), log(summary(:,2)), 'r','LineWidth',2);
        plot(log(summary(:,1)), log(summary(:,4)), '-+b','LineWidth',2);
        %plot(log(summary(:,1)), log(summary(:,3)), '--b','LineWidth',2);
        %plot(log(summary(:,1)), log(summary(:,5)), '--+b','LineWidth',2);
        plot(log(summary(:,1)), log(summary(:,6)), 'black-o', 'LineWidth',2);
        plot(log(summary(:,1)), log(summary(:,7)), ':gs', 'LineWidth',2);

        xlabel('d');
        ylabel('log(mean runtime [s])');
        xticks(log(summary(:,1)))
        xticklabels({'100','200','500','1000','2000','5000','10000'})
        %legend('ROIPCA - 1', 'fROIPCA - 1', 'ROIPCA - 2', 'fROIPCA - 2','IPCA', 'CCIPCA');
        legend('ROIPCA', 'fROIPCA','IPCA', 'CCIPCA');
        
        