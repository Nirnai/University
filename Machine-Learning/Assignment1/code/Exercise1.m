function [par, p1, p2] = Exercise1(k)
%This function computes a polynomial model for a robot control. The control
%inputs are the velocity v and angular velocity w. These map on an x,y
%location and an orientation theta. 
%   Using regression the ideal model parameters are learned. By doing a
%   k-fold crossvalidation the optimal model complexity is determined. The
%   input to the function is the fold number k. The output are the
%   modelparamters A and the complexity p1 for location and p2 for
%   orientation. 

    %%
    load('Data.mat');
    
    % split output data into localization and orientation
    Output1 = Output(1:2,:);
    Output2 = Output(3,:);
    
    [p1,A1] = c_validation(Input, Output1,k);
    [p2,A2] = c_validation(Input, Output2,k);
    
    A = zeros(size(A1,1),3);
    A(:,1:2) = A1;
    A(1:size(A2,1),3) = A2;
    
    latex(sym(vpa(round(A,4))))
    
    par = cell(1,3);
    par{1,1} = A1(:,1);
    par{1,2} = A1(:,2);
    par{1,3} = A2;
    
    save('params','par');
end

%% 

function X = InputMatrix(Input, p)
% This function convertes the input matrix into the form, where the
% variables correspong to the column dimension and the samples to the row
% dimension. Each row of the matrix is composed of the form: 
% row = (1,v,w,vw,v^2,w^2...) 
    X = zeros(length(Input), 4+((p-1)*3));
    X(:,1) = ones(length(Input),1);
    X(:,2:3) = Input';
    X(:,4) = (Input(1,:) .* Input(2,:))';
    
    deg = 1;
    
    while(deg<p)
        X(:,(2+3*deg)) = (Input(1,:).^(deg+1))';
        X(:,(3+3*deg)) = (Input(2,:).^(deg+1))';
        X(:,(4+3*deg)) = ((Input(1,:).*Input(2,:)).^(deg+1))';
        deg = deg + 1;
    end
end

%%
function [A] = polynomial_regression(Input, Output, p)
% Computes a polynomial fit for input-, output-data pairs for 
% a robot navigation application.
%   The Input is given by velocity v and angular velocity w of the robot.
%   The Output is given by the x,y position or theta orientation of the
%   robot. This function finds a plolynomial fit for the data. 
%   The Output of this Function is a Matrix A with the
%   coefficients of the polynomial. The position and orientation each have
%   a own optimal polinomial fit of degree p1 and p2.

    % It is possible to model each of the two polinomials as a matrix
    % multiplication of the kind Ax = y. The number of rows in x depend on
    % the degree of the polynomial. 
    
    X = InputMatrix(Input,p);
    
    Y = Output';
    
    A = pinv(X)*Y;  
end

%%
function[p, A] = c_validation(Input, Output, k)
% This function runs a k-fold crossvalidation on the polynomial regressian
% algorithm implemented earlier. It runs through the polynomials 1 to 6 to
% avoid overfitting and finding just the right degree for the polinomial. 

    data_size = size(Input,2);
    fold_size = floor(data_size/k);
    total_index = 1:data_size;
    fold_index = total_index(1:fold_size:end)-1;
    %fold_idx(end) = data_size;
    
    
    val_error = zeros(1,6);
    for p = 1:6
        error = zeros(1,k);
        for n = 1:k
            % Training Data
            train_input = Input;
            train_input(:,fold_index(n)+1:fold_index(n)+fold_size) = [];
            train_output = Output;
            train_output(:,fold_index(n)+1:fold_index(n)+fold_size) = [];
            
            % Test Data
             test_input = Input(:,fold_index(n)+1:fold_index(n)+fold_size);
             test_output = Output(:,fold_index(n)+1:fold_index(n)+fold_size)';
             
            % predicted parameters
            A = polynomial_regression(train_input, train_output, p);
            % predicted output
            prediction = InputMatrix(test_input,p) * A;
            
            if size(prediction,2) == 2
               error(n) = mean(sqrt((test_output(:,1) - prediction(:,1)).^2+(test_output(:,2) - prediction(:,2)).^2)); 
            elseif size(prediction,2) == 1
                error(n) = mean(sqrt((test_output(:,1) - prediction(:,1)).^2));
            end
        end
        val_error(p) = mean(error);
    end
    
    [e,p] = min(val_error);

    A = polynomial_regression(Input, Output, p);
end




