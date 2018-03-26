A = importdata('A.txt',' ');
B = importdata('B.txt',' ');
p = importdata('pi.txt',' ');
Test = importdata('Test.txt',' ');

%% Forward procedure
% initialization
classif = zeros(1,10);
loglikelihood = zeros(1,10);
for r = 1:10
    o = Test(:,r);
    alpha = zeros(60,12);
    for n = 1:12
        alpha(1,n) = p(1,n)*B(o(1,1),n);
    end

    for i=2:size(o,1)
        for n = 1:12
            alpha(i,n) = sum(alpha(i-1,:) .* A(:,n)',2) * B(o(i,1),n);
        end
    end

    loglikelihood(1,r) = log(sum(alpha(end,:),2));

end
