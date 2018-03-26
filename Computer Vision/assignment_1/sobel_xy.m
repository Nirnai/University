function [Fx,Fy] = sobel_xy(Image)
% In dieser Funktion soll das Sobel-Filter implementiert werden, welches
% ein Graustufenbild einliest und den Bildgradienten in x- sowie in
% y-Richtung zurückgibt.

% horizontal sobel filter Sx
Sx = [1 0 -1; 2 0 -2; 1 0 -1];
% vertical sobel filter Sy
Sy = [1 2 1; 0 0 0; -1 -2 -1];

% perform 2D matrix convolution: F = S * A
Fx = conv2(double(Image), double(Sx), 'same');
Fy = conv2(double(Image), double(Sy), 'same');

end

