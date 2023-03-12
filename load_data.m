function [X, d_, Q_gt, sigma] = load_data(data_type, n, d, m)
    Q_gt = [];
    sigma = 0;
    if strcmp(data_type, 'mnist') == 1
        P = randperm(60000);
        X = loadMNISTImages('Datasets/train-images.idx3-ubyte');X = X';X = X(P,:);
    elseif strcmp(data_type, 'poker') == 1
        X = csvread('Datasets/poker-hand-training-true.data');   X = X(randperm(size(X,1)), 1:(end - 1));
    elseif strcmp(data_type, 'wine') == 1
        X = load('Datasets\wine_dataset.mat'); X = X.wine; X = X(randperm(size(X,1)),:);
    elseif strcmp(data_type, 'har') == 1
        X = load('Datasets\har_x_train.mat'); X = X.harXtrain; X = X(randperm(size(X,1)),:);
    elseif strcmp(data_type, 'super') == 1
        X = load('Datasets\superconductor_dataset.mat'); X = table2array(X.X); X = X(randperm(size(X,1)),:);
    elseif strcmp(data_type, 'gauss_low_rank') == 1
        sigma = 15;
        Q = orth(randn(d, m ));
        D = diag(linspace(20,5,m));
        cov_matrix = Q * D * Q';
        cov_matrix = cov_matrix + sigma * eye(d);
        Q_gt = Q;
        X = mvnrnd(zeros(1,d),cov_matrix, n);
        
        %[~,~,Q_gt] = svds(X, m);
    elseif strcmp(data_type, 'not_low_rank') == 1
        Q = orth(rand(d, m));
        %D = diag(sort(rand(m,1), 'descend'));
        D = diag(linspace(0.6,0.5, m));
        cov_matrix = Q * D * Q' + eye(d);
        
        X = mvnrnd(zeros(1,d),cov_matrix, n);
    elseif strcmp(data_type, 'comp_paper_example') == 1
        for i=1:d
            for j = 1:d
            	cov_matrix(i,j) = min(i,j)/d;
            end
        end
        X = mvnrnd(zeros(1,d),cov_matrix, n);
    elseif strcmp(data_type, 'mult_noise') == 1
        Q = orth(randn(d, m));
        D = diag(linspace(11,7,m));
        cov_matrix = Q * D * Q';
        
        cov_matrix = cov_matrix;
        Q_gt = Q;
        X = mvnrnd(zeros(1,d),cov_matrix, n);
        X = X .* (randn(size(X)) * 2 + 1);
    elseif strcmp(data_type, 'missing_random') == 1
        Q = orth(randn(d, m));
        D = diag(linspace(20,5,m));
        cov_matrix = Q * D * Q';
        sigma = 2;
        cov_matrix = cov_matrix;
        Q_gt = Q;
        X = mvnrnd(zeros(1,d),cov_matrix, n);
        %[~,~,Q_gt] = svds(X,m);
        S = binornd(1,0.5,size(X));
        X = X .* S + randn(size(X)) * sigma;
        %[~,~,Q_gt] = svds(X,m);
    elseif strcmp(data_type, 'outliers_random') == 1
        Q = orth(randn(d, m));
        D = diag(linspace(20,5,m));
        cov_matrix = Q * D * Q';
        
        cov_matrix = cov_matrix;
        Q_gt = Q;
        
        X = mvnrnd(zeros(1,d),cov_matrix, n);
        the_mean = mean(mean(X));
        %[~,~,Q_gt] = svds(X,m);
        S = binornd(1,0.95,size(X)) .* (rand(size(X)) * 100 * the_mean - 50 * the_mean);
        X = X + S;
    elseif strcmp(data_type, 'noise_outliers') == 1
        Q = orth(randn(d, m));
        D = diag(linspace(20,5,m));
        cov_matrix = Q * D * Q';
        
        cov_matrix = cov_matrix;
        Q_gt = Q;
        pct_ouliers = 0.01;
        n_outliers = floor(n * pct_ouliers);
        X = mvnrnd(zeros(1,d),cov_matrix, n);
        P = randperm(n); P = P(1:n_outliers);
        X(P, :) = 0;
        X = X + randn(size(X)) * 3;
        [~,~,Q_gt] = svds(X,m);
    else
        X = zeros(n,d);

    end
    
    X = X(1:n, :);
    %X = zscore(X);
    X = bsxfun(@minus,X,mean(X));
    d_ = size(X,2);
end

