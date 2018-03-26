function  Merkmale = harris_detektor(Image,varargin) 
% In dieser Funktion soll der Harris-Detektor implementiert werden, der
% Merkmalspunkte aus dem Bild extrahiert
% Um die die Funktion der jeweiligen Anwendung anzupassen, k�nnen folgende
% Optionale Parameter �bergeben werden:

    %% Input parser
    P = inputParser;
    
    default_do_plot = false;
    default_segment_length = 3;
    default_k = 0.05;
    default_tau = 1E8;
    default_min_dist = 10;
    % tile_size muss so gew�hlt werden, dass eine ganzzahlige
    % Kachelanzahl in das zu analysierende Bild passen!!!!
    default_tile_size = [100, 100];
    default_N = 10;
    default_sigma = 3;

    % Liste der notwendigen Parameter
    % Ein Bild als Input ist zwingend notwendig
    P.addRequired('Image', @(x)isnumeric(x));

    % Liste der optionalen Parameter
    
    P.addOptional('do_plot', default_do_plot, @(x) islogical(x));
    P.addOptional('segment_length', default_segment_length,  @(x) x >= 3 && mod(x,2));
    P.addOptional('k', default_k, @(x)isnumeric(x));
    P.addOptional('tau', default_tau, @(x)isnumeric(x));
    P.addOptional('min_dist', default_min_dist);
    P.addOptional('tile_size', default_tile_size);
    P.addOptional('N', default_N, @(x)isnumeric(x));
    P.addOptional('sigma', default_sigma, @(x)isnumeric(x));
    
    % Lese den Input
    P.parse(Image,varargin{:});
    
    % Extrahiere die Variablen aus dem Input-Parser
    do_plot = P.Results.do_plot;
    segment_length = P.Results.segment_length;
    k = P.Results.k;
    tau = P.Results.tau;
    min_dist = P.Results.min_dist;
    tile_size = P.Results.tile_size;
    
    if numel(tile_size) == 1
        tile_size   = [tile_size,tile_size];
    end
    
    N = P.Results.N;   
    sigma   = P.Results.sigma;

    
    
    
    %% Filter
    
    % Wenn das Bild nicht in graustufen vorliegt, kann keine analyse
    % stattfinden. �berpr�fung der anzhal der Kan�le
    if size(Image,3) == 1
        
        % Bildgradient in x- und y-Richtung
        [Ix, Iy] = sobel_xy(Image);

        % Berechnung der eintr�ge der Hassis-Matrix
        Ix2 = Ix .* Ix;
        Iy2 = Iy .* Iy;
        Ixy = Ix .* Iy;


        % Gaussian Filter H mit variabler Fenstergr��e
        hfsize = floor(segment_length/2);            % Halbe Fesntergroesse
        ind = -hfsize:hfsize;               % indices (Bsp.: -1,0,1)
        h = exp(-(ind.^2)/(2*sigma^2))';    % 1D Gauss
        H = h*h';                           % 2 Separable Filter ergeben das 2D Filter
        H = H./sum(sum(H));                 % Summe aller Koeffizienten ergibt 1

        % Anwendung des Gau�filters auf Harris Operatoren
        Ix2 = conv2(Ix2, H, 'same');
        Iy2 = conv2(Iy2, H, 'same');
        Ixy = conv2(Ixy, H, 'same');

        % Kriterium f�r Ecken oder Kanten: H = det(G) - k * tr^2(G)
        % Mit diesem vereinfachten Kriterium k�nnen Kanten und Ecken sehr
        % sparsam detektiert werden. Es is allerdings nicht m�glich die
        % orientierung der Kanten zu bestimmen. 
        H = ((Ix2 .* Iy2) - (Ixy .^2)) - k * (Ix2 + Iy2).^2;

        % �berpr�fung auf Eckken. Eine Ecke liegt vor wenn die Harris
        % Matrix einen gro�en Wert annimmt. Grenzwert kann �ber tau und k
        % eingestellt werden.
        H(H < tau) = 0;
        [row,col] = find(H);
        
        
        
        %% Erste Option Merkmale zu extrahieren.
        %Merkmale = [col, row]';
        
        
        
        %% Merkmalsreduktion durch einf�hrung einer mindest Distanz zwischen Merkmalen
        % Erzeugung einer Maske die Elemente innerhalb eines Kreises des gefundenen Merkmals 
        % mit einem radius = min_dist eliminiert. Maske ist eine (2*min_dist)^2 Matrix. 
        mesh_size = (2*min_dist)-1;
        [xdim, ydim] = meshgrid(1:mesh_size);
        center = min_dist;
        mask = ~(sqrt((xdim - center).^2 + (ydim - center).^2) <= min_dist);
        mask(center, center) = 1;
        
        % Um diese Maske auch am Rand anwenden zu k�nnen, muss die Matrix mit
        % min_dist 'gepadded' werden. Dies sorgt f�r eine Indexverschiebung.
        row = row + min_dist;
        col = col + min_dist;
        H = padarray(H,[min_dist,min_dist]);
        
        % Damit daf�r gesorgt ist, dass das st�rkste merkmal die
        % schw�cheren in der n�he eleminiert, und nicht vis versa,
        % werden die Merkmale nach ihrer St�rke sortiert
        lin_idx = sub2ind(size(H), row, col);
        [sorted, true_idx] = sort(H(lin_idx), 'descend');
        % Indizes um Urspr�ngliche Reihenfolge Rekonstruieren zu koennen.
        row = row(true_idx);
        col = col(true_idx);
        
        % Wenn das Merkmal noch nicht elimiert wurde, wird die Maske
        % angewendet um nachbarn zu elimieren, bis der Mindestabstand bei
        % allen Mermalen eerf�llt ist.
        for n = 1:length(sorted)
            if H(row(n), col(n))>0
               row_interval = (row(n) - (min_dist-1)):(row(n)+(min_dist-1));
               col_interval = (col(n) - (min_dist-1)):(col(n)+(min_dist-1));
               H(row_interval, col_interval) = H(row_interval, col_interval) .* mask;
            end
        end

        %% Zweite Option Merkmale zu extrahieren. Es geht jedoch noch besser!
        %[row, col] = find(H);
        %Merkmale = [col - min_dist,row - min_dist]';
       
        %% Merkmalsreduktion durch Einf�hrung von Kacheln und einer Mindestanzahl an Merkmalen pro Kachel zu gew�hrleisten.
        % Vorheriges Padding muss R�ckg�ngig gemacht werden, um wieder die
        % Originalgr��e des Bildes zu haben. 
        H = H(min_dist:end-(min_dist+1), min_dist:end-(min_dist+1));
        vertical =  ones(1,size(H,1)/tile_size(1))*tile_size(1);
        horizontal = ones(1,size(H,2)/tile_size(2))*tile_size(2);
        
        % Einteilung der H Matrix in von 'tile_size' abh�ngig gro�en
        % Kacheln.
        C = mat2cell(H, vertical, horizontal);
        
        % Suchen der N maximalen Werte in einer Kachel und elimieren aller
        % anderen. 
        for n = 1:size(vertical,2)
            for m = 1:size(horizontal,2)
              temp = C{n,m};
              temp = reshape(temp, 1, numel(temp));
              temp = sort(temp,'descend');
              lowest_max = temp(N);
              C{n,m}(C{n,m}<lowest_max)=0;
            end
        end
        % Rekonstruktion der H Matrix mit den modifizierten Kacheln
        H = cell2mat(C);
        
        %% Randbehandlung
        % Die Faltung mit mit dem Gau�filter hat daf�r gesorgt, dass die
        % R�nder nicht ordnungsgem�� in die Faltung mit einbezogen wurden.
        % Dies kann zur fehlerhaften Merkmalserkennung am Rand f�hren. 
        % Da in dieser Implementieren davon ausgegangen wird, dass die
        % Randpixel keine f�r den Anwender interessanten Merkmale
        % aufwei�en, werden diese einfach vernachl�sstig, bzw. alle
        % gefunden Merkmale gel�scht. Dies spart rechenzeit und resourcen. 
        % Falls diese Annahme nicht mehr getroffen werden kann, muss der
        % Rand gesondert behandelt werden. z.B. Um einen sanften �bergang zu
        % gew�hrleisten, k�nnte das zu analysierende Bild an allen R�ndern
        % dupliziertwerden um nach der Faltung wieder entwert zu werden. 
        H(:, 1:segment_length) = 0;
        H(1:segment_length, :) = 0;
        H(:, end-segment_length:end) = 0;
        H(end-segment_length : end, :) = 0;
        
        %% Dritte Option die Merkmale zu extrahieren. Ergebnis ist zufriedenstellen.
        [row, col] = find(H);
        Merkmale = [col,row]';

        
        %% Verschanschaulichung der gefundenen Merkmale im Original Bild.
        if (do_plot) 
            imshow(Image);
            hold on
            plot(Merkmale(1, :), Merkmale(2,:), 'r.');
        end

    else
        disp('Image is not in grayscale format...');
        return;
    end

end