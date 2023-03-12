%{
This script perfroms comprison between various online PCA methods.
It reproduces Table 2 and Table 3 in the paper ROIPCA: An online memory-restricted PCA algorithm based on rank-one
updates

parameters: 
data_type = string of the dataset to use. Options: 'mnist', 'poker',
'wine', super' (for Table 3), 'comp_paper_example' (for Table 2),
'not_low_rank' (for table 4)

n_trials = # of experiemts to perform before taking the median result
n_updates = # of data samples to process in the comparison
n_train = #  of data samples to perform CV in the methods that require it
n0 = size of the initial dataset
d = for synthetic datasets, the diemension of the data
m = # of components to calculate
auto_m = if 1, sets m to account 80% of the varince, over-ridding the m
above
%}

clear;
addpath('../Datasets')  
data_type = 'comp_paper_example';
n_trials = 5;
n_updates = 10000;
n_train = 1000;
n0 = 500;
d = 100;


n = n0 + n_updates;
m = 5;
auto_m = 1;
downsample_rate = 100;

for trial = 1:n_trials

    fprintf('[+] trial: %d\n', trial);
    
    graph_idx = 2;

    [X, d] = load_data(data_type, n0 + n_updates, d, m);
    [Xtrain, ~] = load_data(data_type, n_train, d, m);
    if auto_m
        the_eigs = sort(real(eig(X' * X)), 'descend');
        m = min(10, find(cumsum(the_eigs)/sum(the_eigs) > 0.8, 1));
        fprintf('automatic m = %d\n', m);
    end 
    
    % init parameters
    if 1 
        eta_hebbian = find_eta('Hebbian', Xtrain, Xtrain, [], [], m);
        eta_oja = find_eta('Oja', Xtrain, Xtrain, [], [], m);
        %sizes_dbpca = find_eta('DBPCA', Xtrain, Xtrain, [], [], m);
    end
    
    m_for_error = m;

    % CALC EIGENSPACES
    [Q_final, D_finals] = eigs(X' * X , m);
    Q_final_for_error = Q_final(:, 1:m_for_error);
    norm_ = norm(Q_final_for_error * Q_final_for_error', 'fro')^2; 
    
    X0 = X(1:n0,:);
    [Q_0, D_0] = eigs(X0' * X0, m);
    X = X(n0 + 1: end, :);

    % init
    the_trace =  trace(X0' * X0);
    Q_IPCA = Q_0; D_IPCA = D_0;
    Q_CCIPCA = Q_0; D_CCIPCA = D_0;
    Q_DBPCA = Q_0;
    Q_ROI = Q_0; D_ROI = D_0;
    Q_fast = Q_0; D_fast = D_0;
    Q_hebbian = Q_0;
    Q_oja = Q_0;
       
    for i = 1:n_updates
        update = X(i,:)';
        
        [Q_ROI,D_ROI, the_time] = update_eigenspectrum([],norm(update)^2,update/norm(update),Q_ROI,diag(D_ROI), '112', the_trace, i); timing(i, 1) = the_time;
        
        [Q_fast,D_fast, the_time, n_itr, lr_] = update_eigenspectrum_fast([],norm(update)^2,update/norm(update),Q_fast,diag(D_fast), '112', the_trace, i); timing(i, 2) = the_time;
        
        tic; [Q_IPCA, D_IPCA] = IPCA(Q_IPCA, D_IPCA, update); timing(i, 3) = toc;

        tic; [Q_CCIPCA, D_CCIPCA] = CCIPCA(i + n0 , Q_CCIPCA, D_CCIPCA, update); timing(i, 4) = toc;

        tic; [Q_hebbian] = hebbian(Q_hebbian, update, i + n0, eta_hebbian); timing(i, 5) = toc;
        
        tic; [Q_oja] = OJA(Q_oja, update, i + n0, eta_oja); timing(i, 6) = toc;

        %tic; [Q_DBPCA] = DBPCA(i, sizes_dbpca, Q_DBPCA, X); timing(i, 7) = toc;
        
        
        the_trace = the_trace + norm(update)^2;
        
        if i == 1
            err = norm(Q_final * Q_final' - Q_0*Q_0', 'fro')^2/norm_;
            err_vecs(1, 1) = err;
            err_vecs(1, 2) = err;
            err_vecs(1, 3) = err;
            err_vecs(1, 4) = err;
            err_vecs(1, 5) = err;
            err_vecs(1, 6) = err;
            err_vecs(1, 7) = err;
        else
            if (mod(i,downsample_rate) == 0)
                fprintf('[+] update: %d\n', i);
                err_vecs(graph_idx, 1) = norm(Q_final * Q_final' - Q_ROI*Q_ROI', 'fro')^2/norm_;
                err_vecs(graph_idx, 2) = norm(Q_final * Q_final' - Q_fast*Q_fast', 'fro')^2/norm_;

                
                err_vecs(graph_idx, 3) = norm(Q_final_for_error * Q_final_for_error' - Q_IPCA(:,1:m_for_error)*Q_IPCA(:,1:m_for_error)', 'fro')^2/norm_;
                err_vecs(graph_idx, 4) = norm(Q_final_for_error * Q_final_for_error' - Q_CCIPCA(:,1:m_for_error)*Q_CCIPCA(:,1:m_for_error)', 'fro')^2/norm_;
                err_vecs(graph_idx, 5) = norm(Q_final_for_error * Q_final_for_error' - Q_hebbian(:,1:m_for_error)*Q_hebbian(:,1:m_for_error)', 'fro')^2/norm_;
                err_vecs(graph_idx, 6) = norm(Q_final_for_error * Q_final_for_error' - Q_oja(:,1:min(m_for_error, size(Q_oja,2)))*Q_oja(:,1:min(m_for_error, size(Q_oja,2)))', 'fro')^2/norm_;
                err_vecs(graph_idx, 7) = norm(Q_final_for_error * Q_final_for_error' - Q_DBPCA(:,1:m_for_error)*Q_DBPCA(:,1:m_for_error)', 'fro')^2/norm_;
                graph_idx = graph_idx + 1;
            end
        end
        
    end
    
    err_roipca1_total(trial,:) = err_vecs(:, 1)';
    err_froipca1_total(trial,:) = err_vecs(:, 2)';
    err_ipca_total(trial,:) = err_vecs(:, 3)';
    err_ccipca_total(trial,:) = err_vecs(:, 4)';
    err_hebbian_total(trial,:) =  err_vecs(:, 5)';
    err_oja_total(trial,:) =  err_vecs(:, 6)';
    err_dbpca_total(trial,:) =  err_vecs(:, 7)';
end

% print table of error
fprintf('no update: %.8d+-%.8d \n', median(err_roipca1_total(:, 1)), std(err_roipca1_total(:, 1)));
fprintf('ROIPCA: %.8d+-%.8d \n', median(err_roipca1_total(:, end)), std(err_roipca1_total(:, end)));
fprintf('fROIPCA: %.8d+-%.8d \n', median(err_froipca1_total(:, end)), std(err_froipca1_total(:, end)));
fprintf('IPCA: %.8d+-%.8d \n', median(err_ipca_total(:, end)), std(err_ipca_total(:, end)));
fprintf('CCIPCA: %.8d+-%.8d \n', median(err_ccipca_total(:, end)), std(err_ccipca_total(:, end)));
fprintf('HEBBIAN: %.8d+-%.8d \n', median(err_hebbian_total(:, end)), std(err_hebbian_total(:, end)));
fprintf('OJA: %.8d+-%.8d \n', median(err_oja_total(:, end)), std(err_oja_total(:, end)));
%fprintf('DBPCA: %.8d+-%.8d \n', median(err_dbpca_total(:, end)), std(err_dbpca_total(:, end)));

% plot graph of error
hold on;
len = size(err_roipca1_total,2);
x_axis = downsample_rate * (1:len);
plot(x_axis, median(err_roipca1_total), 'r','LineWidth',1.5);
plot(x_axis, median(err_froipca1_total), '--r', 'LineWidth',1.5);       
plot(x_axis, median(err_ipca_total), 'black', 'LineWidth',1.5);
plot(x_axis, median(err_ccipca_total), 'g', 'LineWidth',1.5);
plot(x_axis, median(err_hebbian_total), 'c', 'LineWidth',1.5);
plot(x_axis, median(err_oja_total), 'm', 'LineWidth',1.5);
%plot(x_axis, median(err_dbpca_total), '-+m', 'LineWidth',1.5);
hold off;
xlabel('# iteration');
ylabel('error');
xlim([0,n_updates]);
legend('ROIPCA', 'fROIPCA', 'IPCA', 'CCIPCA', 'Hebbian', 'Oja');
