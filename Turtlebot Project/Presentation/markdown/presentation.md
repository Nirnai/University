# Control and Navigation
* Bestimmen der Pose aus Odometry

    - x und y direkt aus odom und theta aus quaternion

    - Erweiterung mit gmapping tfs map and base_footprint um drift freie lokalisierung zu bekommen

* Bekommt ein Ziel als x und y coordinate und bestimmt nötige rotation selbständig

    - Entscheidung in welche Richtung gedreht werden soll

* P-Regler für die letzte Meter, um Genauigkeit zu erhöhen

    - Eventuell noch auf PID erweitern

* Erkennt Objekte aus Lidar und stoppt wenn zu nahe

    - Noch zu machen: avoidance path

# Object Detection 
* Sucht kontinuierlich im LIDAR nach Objekten und published deren Posen.

* Führt die Initialisierung durch:
 
    + Es wird nach 3er Objektreihen gesucht, die 2 Bedingungen erfüllen:

        - Das Verhältnis zwischen den 2 inneren Abstanden        beträgt ca. 33%

        - Ale 3 Objekte liegen auf einer Geraden.

    + Berechnung von a/b mittels Kosinussatz

    + Berechnung von der eigenen Position im Spielfeld, Spielfeld und gegnerischem Tor
* Initiale Messung wird etwa 10x durchgeführt
* Je nach Abweichung der Ergebnisse werden sie akzeptiert oder nicht.
* Wenn nicht, alterative Initialisierung durch Bewegung.

# Image Processing 

* Läuft auf dem Turtlebot

* Sucht nach den blau/gelben Pucks und den grünen Pfosten

    - Thresholding der Farben im Bild in HSV

    - Canny Kantendetektor über Threshold Img

    - Konturen Suche mit findContours

    - Berechnung der HuMoments für Konturen ab einer gewissen Größe

    - Bestimmung des Objekts anhand Farbe und der 7 HuMoments

    - Bestimmung des Mittelpunktes der Objekte im Bild

    - Bestimmung der 3D Pose des Objekts aus der Punktwolke

* Bestimmt Teamfarbe

    - Nach der Initialisierung fahren wir seitlich vor das Tor und nehmen ein Bild auf.

    - Im Bild wird nach einer gelben Kontur gesucht und deren Fläche wird approximiert.

# Angelina

* Umwandlung der testgui als in ausführbaren node

* ergänzen eines topics zur Kommunikation mit dem turtlebot