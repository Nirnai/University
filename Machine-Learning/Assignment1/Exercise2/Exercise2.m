function [error_opt, d_opt] = Exercise2( dmax )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    % load data
    img_train = loadMNISTImages('train-images.idx3-ubyte')';
    labels_train = loadMNISTLabels('train-labels.idx1-ubyte');
    img_test = loadMNISTImages('t10k-images.idx3-ubyte')';
    labels_test = loadMNISTLabels('t10k-labels.idx1-ubyte');

    % Preprocess Data
    X_train = img_train - mean(img_train);
    Y_train = labels_train;
    X_test = img_test - mean(img_train);
    Y_test = labels_test;
    
    % compute classification error
    total_clasif = size(Y_test,1);
    error_rate = zeros(1,dmax);
    prediction = zeros(size(X_test,1),60);
    for d = 1:dmax
        fprintf('d = %i\n',d);
        % training
        [X_d,U_d] = train(X_train,d);
        % prediction
        prediction(:,d) = predict(X_test(:,:), X_d, Y_train,U_d);
        diff = prediction(:,d) - Y_test;
        diff(diff~=0) = 1;
        false_clasif = sum(diff);  
        error_rate(1,d) = false_clasif/total_clasif;    
    end
    
    [error_opt,d_opt] = min(error_rate);
    
    % Error plot
    figure()
    title('Classification Error')
    plot(error_rate)
    hold on;
    plot(d_opt, error_opt,'o','MarkerSize',10);
    xlabel('d Principal Components');
    ylabel('error rate');
    
    % Confusion Matrix
    conf_mat = confusionmat(Y_test, prediction(:,d_opt));
    helperDisplayConfusionMatrix(conf_mat);
    
end


function [X_d,U_d] = train(X,d)
    C = cov(X);
    [U,D] = eig(C);
    [D,I] = sort(diag(D),'descend');
    U_d = U(:,I);
    U_d = U_d(:,1:d);
    X_d = X * U_d;
end

function [prediction] = predict(test,train_d,train_labels,U_d)
    X_d = test * U_d;
    likelihood = zeros(size(test,1),10);
    for class=0:9
        [mu, Sigma] = gaussian_parameters(train_d, train_labels,class);
        likelihood(:,class+1) = mvnpdf(X_d, mu, Sigma);
    end
    [M,I] = sort(likelihood,2,'descend');
    I = I(:,1); 
    prediction = I - 1;
end

function [mu,Sigma] = gaussian_parameters(Input,Labels,Class)
    In = Input(Labels == Class, :);
    mu = mean(In);
    Sigma = cov(In);
end

