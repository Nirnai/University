function [] = Exercise3_kmeans()
% This function loads the gesture dataset provided in the folder and
% computes 7 cluster centers for each gesture.
%   It does so by using the k-means algorithm and the given initialization
%   centroids.

load('gesture_dataset.mat');

data1 = reshape(gesture_l,600,3);
data2 = reshape(gesture_o,600,3);
data3 = reshape(gesture_x,600,3);

centroids1 = init_cluster_l;
centroids2 = init_cluster_o;
centroids3 = init_cluster_x;

[labels1, centroids1] = k_means(data1, centroids1, 10e-6);
plot_clusters(data1, labels1, centroids1,1);

[labels2, centroids2] = k_means(data2, centroids2, 10e-6);
plot_clusters(data2, labels2, centroids2,2);

[labels3, centroids3] = k_means(data3, centroids3, 10e-6);
plot_clusters(data3, labels3, centroids3,3);


end

function[] = plot_clusters(data, labels, centroids,plot_nbr)
% This function plots every datapoint belonging to the same cluster center
% in the same color. It also plots the centroids as a filled black circle.
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

function[labels, centroids] = k_means(data, centroids, decrement_threshold)
% This function takes data and initial centroids, and computes new
% centroids by using k_means algorithm. The output are the new- centroids
% and labels
    total_distortion = 1;
    while(1)
        labels = closest_centroid(data, centroids);
        centroids = new_centroids(data, labels);
        new_destortion = evaluate(data, labels, centroids);
        if(abs(new_destortion - total_distortion) <= decrement_threshold)
            break;
        end
        total_distortion = new_destortion;
    end
end

function[labels] = closest_centroid(data, centroids)
    distance = zeros(600,7);
    for k = 1:7
        distance(:,k) = sqrt(sum((data - centroids(k,:)).^2,2));
    end
    
    [m,labels] = min(distance,[],2);
end

function[new_centroids] = new_centroids(data, labels)
    new_centroids = zeros(7,3);
    for k = 1:7
        idx = find(labels == k); 
        new_centroids(k,:) = mean(data(idx,:));
    end
end

function[total_distortion] = evaluate(data,labels,new_centroids)
    J = zeros(1,7);
    for k=1:7
        J(1,k) = sum(sum((data(labels==k) - new_centroids(k,:)).^2,2));
    end
    total_distortion = sum(J);
end