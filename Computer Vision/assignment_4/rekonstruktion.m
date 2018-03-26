function [T,R, lambdas, P1] = rekonstruktion(T1,T2,R1,R2, Korrespondenzen_robust, K)
% Funktion zur Bestimmung der korrekten euklidischen Transformation, der
% Tiefeninformation und der 3D Punkte der Merkmalspunkte in Bild 1

%%
% Extrahiere homogene Merkmalspunkte (Korrespondenzen sind
% Pixelkoordinaten)
if size(Korrespondenzen_robust,1)==4
    x1_pixel = [Korrespondenzen_robust(1:2,:);ones(1,length(Korrespondenzen_robust))];
    x2_pixel = [Korrespondenzen_robust(3:4,:);ones(1,length(Korrespondenzen_robust))];
end

x1 = K^-1 * x1_pixel;
x2 = K^-1 * x2_pixel;

%%
% Erstelle Matrix M
M1 = zeros(1,size(Korrespondenzen_robust,2)+3);
M2 = zeros(1,size(Korrespondenzen_robust,2)+3);
for i=1:size(Korrespondenzen_robust,2)
    
    x2_hat = [0 -1 x2(2,i); 1 0 -x2(1,i); -x2(2,i) x2(1,i) 0];
    
    % Gleichungssystem for R1 und T1
    m1 = zeros(3, size(Korrespondenzen_robust,2)+3);
    m1_1 = x2_hat * R1 * x1(:,i);
    m1_2 = x2_hat * T1;
    m1(:,i) = m1_1;
    m1(:,end-2:end) = m1_2;
    M1 = cat(1,M1,m1);
    
    % Gleichungssystem für R2 und T2
    m2 = zeros(3, size(Korrespondenzen_robust,2)+3);
    m2_1 = x2_hat * R2 * x1(:,i);
    m2_2 = x2_hat * T2;
    m2(:,i) = m2_1;
    m2(:,end-2:end) = m2_2;
    M2 = cat(1,M2,m2);
    
    
end

M1 = M1(2:end,:);
M2 = M2(2:end,:);
%%
[U1,S1,V1] = svd(M1);
[U2,S2,V2] = svd(M2);

lambdas1 = V1(:,end);
lambdas2 = V2(:,end);

test1 = lambdas1;
test1(test1 > 0) = 0;
test1(test1 < 0) = 1;

test2 = lambdas2;
test2(test2 > 0) = 0;
test2(test2 < 0) = 1;


if(sum(test1) > sum(test2))
    lambdas = lambdas2(1:end-3);
    R = R2;
    T = T2;
else
    lambdas = lambdas1(1:end-3);
    R = R1;
    T = T1;
end
P1 = lambdas' .* x1;    

end