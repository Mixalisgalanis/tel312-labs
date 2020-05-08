close all; clear all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Exercise 1 - GRAYSCALE CONVERSION              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Image Resolution is 428 x 320
[X1, map1] = imread('image.png');

%Creating new grayscale map
grayscale_map = zeros;

%Applying Grayscale effect
for i=1:length(map1)
    R = map1(i,1); G = map1(i,2); B = map1(i,3); %Acquiring Color Values
    Y = (222*R + 707*G + 71*B)/1000; %Generating shade of gray
    for j = 1:3
        grayscale_map(i,j) = Y; %Applying shades of gray to grayscale_map
    end
end

%Generates Grayscale Image
imwrite(X1, grayscale_map, 'image_grayscale.png'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Exercise 2 - HSI TO RGB MODEL CONVERSION       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%LAB31239720
k = 7 + 2 + 0;  %k = 9
S = (1-k/10);   %S = 0.1
I = 0.5;

H = zeros;
j = 0;
for i=0:255
    H(i+1) = j;
    j = j + 360/255; 
end

%Determining Max & min values
M = I*(1 + S);
m = 2*I - M;

map_2 = zeros; %256 x 3
%Generating RGB Values for each H
for i=1:length(H)
    if (H(i) < 60)
        R = m + (M - m) * (H(i)/60);
        G = m;
        B = M;
    elseif (H(i) < 120)
        R = m;
        G = m;
        B = m + (M - m)*((120-H(i))/60);
    elseif (H(i) < 180)
        R = m;
        G = m + (M - m) * ((H(i)-120)/60);
        B = m;
    elseif (H(i) < 240)
         R = m + (M - m)*((240-H(i))/60);
         G = M;
         B = m;
    elseif (H(i) < 300)
         R = m;
         G = M;
         B = m + (M - m) * ((H(i)-240)/60);
    elseif (H(i) < 360)
         R = m;
         G = m + (M - m)*((360-H(i))/60);
         B = M;
    end
    
    map_2(i,1) = R;
    map_2(i,2) = G;
    map_2(i,3) = B;
end

%Creating RGB Color Palette
X_new = zeros(256); %256 x 256
for i=1:length(X_new)
    X_new(i,:)=i;
end

%Generating image based on RGB values
imwrite(X_new, map_2, 'paletta.bmp');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Exercise 3 - IMAGE DESATURATION                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Using map from 1st exercise
[X, map] = imread('image2.png');
map_3 = rgb2hsv(map);

while (mean(map_3(:,2)) > 0)
    figure;
    imshow(X, map)
    for i = 1:length(map_3)
        %map_3(i,2) is saturation channel of HSI model
        if (map_3(i,2) > 0.2)
            map_3(i,2) = map_3(i,2) - 0.2;
        else
            map_3(i,2) = 0;
        end
    end
    map = hsv2rgb(map_3);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Exercise 4 - White Balance Technique           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[X_4, map_4] = imread('image.png');

if ~isempty(map_4)
    Im = ind2rgb(X_4,map_4);
end

redChannel = Im(:, :, 1);
greenChannel = Im(:, :, 2);
blueChannel = Im(:, :, 3);

meanR = mean2(redChannel);
meanG = mean2(greenChannel);
meanB = mean2(blueChannel);

redChannel = ((redChannel) * meanG / meanR);
greenChannel = ((greenChannel));
blueChannel = ((blueChannel) * meanG / meanB);

% Recombine separate color channels into a single, true color RGB image.
Ima = cat(3, redChannel, greenChannel, blueChannel);

imwrite(Ima,'image_whitebalance.png')
