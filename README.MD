# Sammlung Externer Toolboxen

Dieses Repo enthält eine Sammlung von externen Toolboxen und Funktionen für Matlab. Diese Funktionen werden häufig gebraucht und stammen von diversen Autoren außerhalb der Universität Hannover.

Moritz Schappler, schappler@imes.uni-hannover.de, 2018-03  
(C) Institut für mechatronische Systeme, Universität Hannover

## Benutzung

Für die Benutzung des Repos muss die Pfadinitialisierung durchgeführt werden mit dem Skript `matlab_ext_path_init.m` aus dem Hauptordner.

### Modultests

Die einzelnen Funktionen des Repos können über die Modultests aus dem Ordner `examples_tests` gestartet werden. Die Funktionsfähigkeit des Gesamt-Repos kann mit dem Gesamt-Testskript `matlab_ext_test_repo.m` aus dem Hauptordner geprüft werden.

## Inhalt

Hier erfolgt eine kurze Übersicht der im Repo enthaltenen Funktionen. Die Funktionsweise ist in den jeweiligen Funktionsdateien und Modultests erklärt.

* `Advanced_Setpoints`: Feedforward Motion Control Toolbox von Paul Lambrechts. Enthält analytische Trapez-Profile
* `export_fig`: Dient zum schönen Exportieren von Matlab Figures in hoher Auflösung und als Vektorgrafik. Wird durch Matlab-interne Routinen nur unzureichend abgedeckt.
* `geom2d`, `geom3d`: Funktionen zum Zeichnen und Berechnen von zwei- und dreidimensionalen geometrischen Objekten
* `rvctools`: Robotics Toolbox von Peter Corke. Enthält allgemeine Funktionen zur Berechnung von Robotern.

## Erweiterung

Weitere nützliche externe Funktionen und Toolboxen können ergänzt werden und damit allen Institutsangehörigen zur Verfügung gestellt werden. Da es sich um externe Programme handelt, ist es nur sinnvoll, wenn diese für viele Projekte eingesetzt werden oder eine Abhängigkeit für andere Matlab-Toolboxen darstellen.
