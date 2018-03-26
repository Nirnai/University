%  Gruppennummer: M10
%  Gruppenmitglieder: Lukas Köpf, Fabian Lickert, Holger Müller, Nirnai Rao, Simon Winkler

%% Hausaufgabe 3
%  Bestimmung von robusten Korrespondenzpunktpaaren mittels RANSAC und
%  Berechnung der Essentiellen Matrix.

%  Für die letztendliche Abgabe bitte die Kommentare in den folgenden Zeilen
%  enfernen und sicherstellen, dass alle optionalen Parameter über den
%  entsprechenden Funktionsaufruf fun('var',value) modifiziert werden können.


%% Bilder laden
Image1 = imread('szeneL.jpg');
IGray1 = rgb_to_gray(Image1);

Image2 = imread('szeneR.jpg');
IGray2 = rgb_to_gray(Image2);

%% Harris-Merkmale berechnen
Merkmale1 = harris_detektor(IGray1,'segment_length',9,'k',0.05,'min_dist',50,'N',20,'do_plot',false);
Merkmale2 = harris_detektor(IGray2,'segment_length',9,'k',0.05,'min_dist',50,'N',20,'do_plot',false);



%% Korrespondenzschätzung
tic;
Korrespondenzen = punkt_korrespondenzen(IGray1,IGray2,Merkmale1,Merkmale2,'min_corr',0.92,'do_plot',false);
zeit_korrespondenzen = toc;
disp(['Es wurden ' num2str(size(Korrespondenzen,2)) ' Korrespondenzpunktpaare in ' num2str(zeit_korrespondenzen) 's gefunden.'])



%% Finde robuste Korrespondenzpunktpaare mit Hilfe des RANSAC-Algorithmus
tic;
Korrespondenzen_robust = F_ransac(Korrespondenzen);
zeit_robuste_korrespondenzen = toc;
disp(['Es wurden ' num2str(size(Korrespondenzen_robust,2)) ' robuste Korrespondenzpunktpaare in ' num2str(zeit_robuste_korrespondenzen) 's mittels RanSaC gefunden.'])
% Zeige die robusten Korrespondenzpunktpaare für Debugging
figure;
imshow(IGray1)
hold on
plot(Korrespondenzen_robust(1,:),Korrespondenzen_robust(2,:),'r*')
imshow(IGray2)
alpha(0.5);
hold on
plot(Korrespondenzen_robust(3,:),Korrespondenzen_robust(4,:),'g*')
for i=1:size(Korrespondenzen_robust,2)
    hold on
    x1 = [Korrespondenzen_robust(1,i), Korrespondenzen_robust(3,i)];
    x2 = [Korrespondenzen_robust(2,i), Korrespondenzen_robust(4,i)];
    line(x1,x2);
end
hold off
% Debugging
%% Berechnung der Essentiellen Matrix
load('K.mat');
E = achtpunktalgorithmus(Korrespondenzen_robust,K)  %extra ohne Semikolon, damit es ausgegeben wird (siehe Aufgabenstellung)

