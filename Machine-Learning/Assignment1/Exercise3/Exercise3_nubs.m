function [] = Exercise3_nubs(K)
% This function loads the gesture dataset provided in the directory and
% computes cluster centers
%   To do so it uses the non uniform binary split (nubs) algorithm to
%   compute the centroids. It takes as input the number of desired
%   clusters K.
load('gesture_dataset.mat');

data1 = reshape(gesture_l,600,3);
data2 = reshape(gesture_o,600,3);
data3 = reshape(gesture_x,600,3);

v = [0.08, 0.05, 0.02];

[y, labels] = nubs(data1, K, v);
plot_clusters(data1, labels, y, 1);

[y, labels] = nubs(data2, K, v);
plot_clusters(data2, labels, y, 2);

[y, labels] = nubs(data3, K, v);
plot_clusters(data3, labels, y, 3);

end



function[] = plot_clusters(data, labels, centroids,plot_nbr)
    colors = 'bkrgmyc';

    for k = 1:7
        figure(plot_nbr);
        hold on;
        new_data = data(labels == k,:);
        scatter3(new_data(:,1),new_data(:,2),new_data(:,3),'.',colors(k));
        hold on;
        scatter3(centroids(:,1), centroids(:,2), centroids(:,3),'O', 'k', 'filled');
    end

end

function [y, labels] = nubs(data, K, v)
    % all data is in class 1 in the beginning
    k = 1;
    labels = ones(600,1);
    y = zeros(k,3);

    while(k<K)
        % colculate centroids
        y = update_codevectors(data, labels, y);
        % determine the worst cluster and split it
        worst_centroid = worst_distortion(data, labels, y);
        % two new cluster centers
        y(worst_centroid,:) = y(worst_centroid,:) + v;
        y = cat(1, y, y(worst_centroid,:) - v);
        % allcote data to new centroids
        labels = new_labels(data, y);
        k = k + 1;
    end
    y = update_codevectors(data, labels, y);
end

function[y_new] = update_codevectors(data, labels, y)
    y_new = zeros(size(y));
    for k = 1:size(y,1)
        y_new(k,:) = mean(data(labels == k,:));
    end
end

function[labels] = new_labels(data, y)
    distance = zeros(size(data,1),size(y,1));
    for k = 1:size(y,1)
        distance(:,k) = sqrt(sum((data - y(k,:)).^2, 2));
    end
    [M,I] = min(distance,[] ,2);
    labels = I;
end 

function [I] = worst_distortion(data, labels, y)
    distortion = zeros(1, size(y,1));
    for k = 1:size(y,1)
        distortion(:,k) = sum(sqrt(sum((data(labels == k) - y(k,1)).^2,2)));
    end
    [M,I] = max(distortion,[],2);
end
