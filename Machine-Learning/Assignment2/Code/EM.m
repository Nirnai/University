data = load('dataGMM.mat');
data = data.Data';


%% Initialization
idx = kmeans(data,4);

w = zeros(4,1);
mu = zeros(4,2);
sigma = zeros(2,2,4);

for i=1:4
    w(i,1) = size(data(idx == i,:),1)/size(data,1);
    mu(i,:) = mean(data(idx == i,:),1);
    sigma(:,:,i) = cov(data(idx == i,:));
end

temp = 1e3;
threshold = 1;
while(threshold > 1e-3)
    %% E step
    % probability of each point for each gaussian
    posterior = zeros(300,4);
    denominator = zeros(300,1);
    for i = 1:4
       posterior(:,i) = w(i,1) .* mvnpdf(data, mu(i,:), sigma(:,:,i));
       denominator = denominator + posterior(:,i);
    end

    posterior = posterior ./ denominator;

    %% M step
    for i = 1:4
        w(i,1) = sum(posterior(:,i))/300;
        mu(i,:) = sum(posterior(:,i) .* data, 1) / sum(posterior(:,i));
        sigma(:,:,i) = (posterior(:,i)' .* (data - mu(i,:))') * (data - mu(i,:)) / (sum(posterior(:,i)));
    end

    %% Eval convergance
    l = zeros(300,4);
    for i = 1:4
       l(:,i) = w(i,1) .* mvnpdf(data, mu(i,:), sigma(:,:,i));
    end 
    loglikelihood = sum(log(sum(l,2)),1);
    threshold = loglikelihood - temp;
    temp = loglikelihood;
    
end




%%
figure()
x1 = -0.08 : 0.001 : 0.1;
x2 = -0.1 : 0.001: 0.1;
[X,Y] = meshgrid(x1, x2);
scatter(data(:,1),data(:,2),2,'k');
hold on; %scatter(mu(1,:),mu(2,:),10);
for i = 1:4
    Z = mvnpdf([X(:) Y(:)],mu(i,:),sigma(:,:,i));
    Z = reshape(Z,size(X));
    contour(X,Y,Z)
    hold on; axis equal 
end


