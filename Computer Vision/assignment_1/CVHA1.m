%  Gruppennummer: M10
%  Gruppenmitglieder: Lukas Köpf, Fabian Lickert, Holger Müller, Nirnai Rao, Simon Winkler

clearvars;

%% Hausaufgabe 1
%  Einlesen und Konvertieren von Bildern sowie Bestimmung von 
%  Merkmalen mittels Harris-Detektor. 

%% Bild laden
Image = imread('szene.jpg');
IGray = rgb_to_gray(Image);

%% Harris-Merkmale berechnen
tic;
M  = harris_detektor(IGray,'do_plot',true);
toc;

fprintf('%i Merkmale wurden gefunden.\n' ,length(M));

        
    
