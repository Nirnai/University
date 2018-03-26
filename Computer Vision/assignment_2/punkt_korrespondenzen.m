function [Korrespondenzen] = punkt_korrespondenzen(I1,I2,Mpt1,Mpt2,varargin)
% In dieser Funktion sollen die extrahierten Merkmalspunkte aus einer
% Stereo-Aufnahme mittels NCC verglichen werden um Korrespondenzpunktpaare
% zu ermitteln.


%% Input parser
P = inputParser;

% Liste der optionalen Parameter
% L?nge des Bildsegments
P.addOptional('window_length', 7, @(x) x >= 3 && mod(x,2))
% Gewichtung f?r das Harris-Kriterium
P.addOptional('min_corr', 0.95, @(x) x > 0 && x <= 1);
% Plot ein/aus
P.addOptional('do_plot', false, @islogical);

% Lese den Input
P.parse(varargin{:});

% Extrahiere die Variablen aus dem Input-Parser
window_length   = P.Results.window_length;
min_corr        = P.Results.min_corr;
do_plot         = P.Results.do_plot;


%% Alles auf double casten
I1           = double(I1);
I2           = double(I2);
Mpt1           = double(Mpt1);
Mpt2           = double(Mpt2);


%%
% Mindestabstand vom Reature zum Rand
edge_check = floor(window_length/2);

% Speicherplatz festlegen die Bildsegmente von Bild 1 und 2
patches1 = zeros(window_length^2,size(Mpt1,2));
patches2 = zeros(window_length^2,size(Mpt2,2));

% extrahieren der Bildpatches und Vektorisieren
for n = 1:size(Mpt1,2)
    if (Mpt1(1, n) > edge_check && Mpt1(1, n) <= (size(I1, 2) - edge_check) && Mpt1(2, n) > edge_check && Mpt1(2, n) <= (size(I1, 1) - edge_check))
        y_idx = Mpt1(2,n)- floor(window_length/2) : Mpt1(2,n)+ floor(window_length/2);
        x_idx = Mpt1(1,n)- floor(window_length/2) : Mpt1(1,n)+ floor(window_length/2);
        patch = reshape(I1(y_idx, x_idx),window_length^2,1);
        patches1(:,n) = patch;
    else
        continue;
    end
end

for n = 1:size(Mpt2,2)
    if(Mpt2(1, n) > edge_check && Mpt2(1, n) <= (size(I2, 2) - edge_check) && Mpt2(2, n) > edge_check && Mpt2(2, n) <= (size(I2, 1) - edge_check))
        y_idx = Mpt2(2,n)- floor(window_length/2) : Mpt2(2,n)+ floor(window_length/2);
        x_idx = Mpt2(1,n)- floor(window_length/2) : Mpt2(1,n)+ floor(window_length/2);
        patch = reshape(I2(y_idx, x_idx),window_length^2,1);
        patches2(:,n) = patch;
    else
        continue;
    end
end

% Normieren der Bild Patches
patches1 = (patches1 - mean(patches1))./std(patches1);
patches2 = (patches2 - mean(patches2))./std(patches2);

% Normalized Cross Correlation der Patches beider Bilder berechnen
NCC = (patches2' * patches1)./(window_length^2 - 1);

% Merkmal extrahieren wenn NCC über dem treshold liegt
[i,j] = find(NCC > min_corr);

% Merkmalskoordinaten der Merkmale speichern
p1 = Mpt1(:,j);
p2 = Mpt2(:,i);

Korrespondenzen = cat(1,p1,p2);


% plot
if do_plot && size(Korrespondenzen, 2) > 0
    imshowpair(I1, I2, 'montage');
    hold on;

    plot(p1(1,:),p1(2,:), 'rs');
    plot(p2(1,:)+3000,p2(2,:), 'gx');
    plot([p1(1,:);p2(1,:)+3000],[p1(2,:);p2(2,:)], '-r');
    axis('off')
end

end

