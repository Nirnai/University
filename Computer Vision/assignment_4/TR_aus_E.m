function [T1,R1,T2,R2] = TR_aus_E(E)
% In dieser Funktion sollen die moeglichen euklidischen Transformationen
% aus der Essentiellen Matrix extrahiert werden

[U,S,V] = svd(E);

if(det(E) < 0)
    I = [1 0 0; 0 1 0; 0 0 -1];
    U = U * I;
end



R_z1 = [0 -1 0; 1 0 0;  0 0 1];
R_z2 = [0 1 0; -1 0 0;  0 0 1];

R1 = U * R_z1 * V';
R2 = U * R_z2 * V';
T1 = U * R_z1 * S * U';
T2 = U * R_z2 * S * U';


end