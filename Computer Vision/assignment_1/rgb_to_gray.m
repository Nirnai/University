function [Gray_image] = rgb_to_gray(Image)
% Diese Funktion soll ein RGB-Bild in ein Graustufenbild umwandeln. Falls
% das Bild bereits in Graustufen vorliegt, soll es direkt zurückgegeben werden.

if size(Image,3) == 1 
    Gray_image = Image;
else
    R = Image(:,:,1);
    G = Image(:,:,2);
    B = Image(:,:,3);

    Gray_image = 0.299 * R + 0.587 * G + 0.114 * B; 
end