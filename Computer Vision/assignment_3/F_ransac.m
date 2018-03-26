function [Korrespondenzen_robust] = F_ransac(Korrespondenzen,varargin)
% Diese Funktion implementiert den RANSAC-Algorithmus zur Bestimmung von
% robusten Korrespondenzpunktpaaren

%% Input parser
P = inputParser;

% Liste der optionalen Parameter
% Fensterlänge
P.addOptional('epsilon', 0.5, @isnumeric)
% Minimal geforderte Korrelation
P.addOptional('p', 0.999, @(x)isnumeric(x) && x > 0 && x < 1);
% Plot ein/aus
P.addOptional('tolerance', 0.075, @isnumeric);

% Lese den Input
P.parse(varargin{:});

% Extrahiere die Variablen aus dem Input-Parser
epsilon     = P.Results.epsilon;
p           = P.Results.p;
tolerance   = P.Results.tolerance;

%% RanSaC

% Berechnung der Iterationszahl s
s = ceil(log(1-p) / log(1-(1-epsilon)^8));


% sampson distanz
% Korrespondenzpunkte a und b zu homogenen Koordinaten x,y,1 erweitern
pt1_hom = ones(3, size(Korrespondenzen, 2));
pt1_hom(1:2,:) = Korrespondenzen(1:2, :);

pt2_hom = ones(3, size(Korrespondenzen, 2));
pt2_hom(1:2,:) = Korrespondenzen(3:4, :);

e3 = zeros(3);
e3(2,1) = 1;
e3(1,2) = -1;

% size of current consensus set
temp = 0;
% best recent parameters
F = zeros(3);

% Führe s Iterationen durch
for i = 1:s
    rand_idx = randsample(1:size(Korrespondenzen,2),8);
    rand_korrespondenzen = Korrespondenzen(:, rand_idx);
    EF = achtpunktalgorithmus(rand_korrespondenzen);
    
    e3 = zeros(3);
    e3(2,1) = 1;
    e3(1,2) = -1;

    d_sampson = (diag(pt2_hom' * EF * pt1_hom).^2)' ./ (sum((e3*EF*pt1_hom).^2,1) + sum((pt2_hom'*EF*e3).^2,2)');
    consensus_set = d_sampson(d_sampson < tolerance);
    
    if(size(consensus_set,2)>temp)
        temp = size(consensus_set,2);
        F = EF;
        Korrespondenzen_robust = Korrespondenzen(:,d_sampson < tolerance);
    end

end

end