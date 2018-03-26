function [EF] = achtpunktalgorithmus(Korrespondenzen,K)
% Diese Funktion berechnet die Essentielle Matrix oder Fundamentalmatrix
% mittels 8-Punkt-Algorithmus, je nachdem, ob die Kalibrierungsmatrix 'K'
% vorliegt oder nicht

% Korrespondenzpunkte a und b zu homogenen Koordinaten x,y,1 erweitern
a = ones(3, size(Korrespondenzen, 2));
a(1:2,:) = Korrespondenzen(1:2, :);

b = ones(3, size(Korrespondenzen, 2));
b(1:2,:) = Korrespondenzen(3:4, :);

% Falls Kalibrierungsmatrix vorhanden, kalibriere Koordinaten mit K
if (nargin > 1)
    % a_kalibriert = K^(-1) * a_unkalibriert
    a = K\a;
    % b_kalibriert = K^(-1) * b_unkalibriert
    b = K\b;
end

% Berechne Kronecker-Produkte und bilde Matrix A
for i = 1:size(Korrespondenzen, 2)
    A(i, :) = transpose(kron(a(1:3, i), b(1:3, i)));
end

% Löse Minimierungsproblem mittels SVD
[U, S, V] = svd(A);

% Bilde G^S aus 9. Spalte von V, G wird aus G^S durch umsortieren gewonnen
G = reshape(V(:, end), [3,3]);

% Führe erneut SVD durch
[Ug, Sg, Vg] = svd(G);

if (nargin > 1) % wenn K gegegeben
    % Kalibrierte Kamera, schätze Essentielle Matrix E
    EF = Ug * [1 0 0; 0 1 0; 0 0 0] * transpose(Vg);
else % wenn kein K gegegeben
    % Unkalibrierte Kamera, setze dritten Singulärwert 0
    Sg(3,3) = 0;
    % berechne Fundamentalmatrix
    EF = Ug * Sg * transpose(Vg);
end

end