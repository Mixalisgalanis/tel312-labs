%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BASIC INFORMATION                                   %
% Course: Digital Image Processing - Lab 2            %
% Deadline: 27-03-2019                                %
% LAB31239720:  Pantelis Karamailis, 2016030040       %
%               Kostantinos Vlachos, 2016030041       %
%               Mixalis Galanis,     2016030036       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Clearing things up
close all;
clearvars;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Exercise 1 - HISTOGRAM, EQUALIZATION, B-W CONVERSION %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1.1 - Displaying Histogram of an overexposed image
Initial_image_1 = imread('overexposed_image.jpg'); %Gathering Image (Resolution: 3960x2640)
% 1.2 - Equalizing image
Equalized_image = histeq(Initial_image_1);    %Applying Equalizer
figure(1);
subplot(2,2,1); imshow(Initial_image_1); title("IMAGE (BEFORE)");   %Displaying Before Image on screen
subplot(2,2,2); imshow(Equalized_image); title("IMAGE (AFTER)");   %Displaying After Image on screen
subplot(2,2,3); imhist(Initial_image_1); title("HISTOGRAM (BEFORE)");	%Displaying Before Histogram on screen
subplot(2,2,4); imhist(Equalized_image); title("HISTOGRAM (AFTER)");   %Displaying After Histogram on screen
% 1.3 - Converting image to Black & White
[counts,binLocations] = imhist(Initial_image_1); %Gathering Histogram Data
%Finding Threshold
half_channels = sum(counts)/2;
sum_of_counts = 0; Threshold = 0;
for i = 1:length(counts)
   sum_of_counts = sum_of_counts + counts(i);
   if (sum_of_counts >= half_channels)
       Threshold = i;
       break;
   end
end
%Conversion
bw_image_1 = rgb2gray(Initial_image_1);
for i = 1 : length(bw_image_1(:,1))
    for j = 1 : length(bw_image_1(1,:))
        if (bw_image_1(i,j) < Threshold)
            bw_image_1(i,j) = 0;
        else
            bw_image_1(i,j) = 255;
        end
    end
end
imwrite(bw_image_1, 'bw_image.jpg');  %Saving image to drive
%Showing side by side comparison
figure(2); 
subplot(2,1,1); imshow(Initial_image_1); title("IMAGE (BEFORE)");   %Displaying Before Image on screen
subplot(2,1,2); imshow(bw_image_1); title("BLACK & WHITE IMAGE (AFTER)");%Displaying After Image on screen

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Exercise 2 - MEDIAN & GAUSSIAN FILTERS               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Initial_image_2 = imread('brain.gif'); %Gathering Image (Resolution: 256x256)
k1 = 7 + 2 + 0;     %k1 = 9
k2 = k1 + 2;        %k2 = 11
%Applying median filters
Filtered_image_2_1 = medfilt2(Initial_image_2, [k1 k1]);
Filtered_image_2_2 = medfilt2(Initial_image_2, [k2 k2]);
%Applying gaussian filters
h_1 = fspecial('gaussian', [k1 k1]);
h_2 = fspecial('gaussian', [k2 k2]);
Filtered_image_2_3 = imfilter(Initial_image_2, h_1);
Filtered_image_2_4 = imfilter(Initial_image_2, h_2);
%Displaying and comparing images
figure(3);
subplot(2,1,1); imshow(Initial_image_2); title("ORIGINAL IMAGE");
subplot(4,2,5); imshow(Filtered_image_2_1); title("IMAGE WITH MEDIAN FILTER [k1 k1]");
subplot(4,2,6); imshow(Filtered_image_2_2); title("IMAGE WITH MEDIAN FILTER [k2 k2]");
subplot(4,2,7); imshow(Filtered_image_2_3); title("IMAGE WITH GAUSSIAN FILTER [k1 k1]");
subplot(4,2,8); imshow(Filtered_image_2_4); title("IMAGE WITH GAUSSIAN FILTER [k2 k2]");


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Exercise 3 - NOISE & NOISE REDUCTION USING FILTERS   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Initial_image_3 = imread('brain.gif'); %Gathering Image (Resolution: 256x256)
%k1, k2 remain the same
%Adding Noise to the Image
Noise_image_3 = imnoise(Initial_image_3,'salt & pepper', 0.3); % very high noise density
%Applying noise reduction filter using median
Filtered_noise_image_3_1 = medfilt2(Noise_image_3, [k1 k1]);
Filtered_noise_image_3_2 = medfilt2(Noise_image_3, [k2 k2]);
%Applying noise reduction filter using average
h_3 = fspecial('average', [k1 k1]);
h_4 = fspecial('average', [k2 k2]);
Filtered_noise_image_3_3 = imfilter(Noise_image_3, h_3);
Filtered_noise_image_3_4 = imfilter(Noise_image_3, h_4);
%Displaying and comparing images
figure(4);
subplot(2,2,1); imshow(Initial_image_3); title("ORIGINAL IMAGE");
subplot(2,2,2); imshow(Noise_image_3); title("IMAGE WITH NOISE");
subplot(4,2,5); imshow(Filtered_noise_image_3_1); title("DENOISED IMAGE WITH MEDIAN FILTER - [k1 k1]");
subplot(4,2,6); imshow(Filtered_noise_image_3_2); title("DENOISED IMAGE WITH MEDIAN FILTER - [k2 k2]");
subplot(4,2,7); imshow(Filtered_noise_image_3_3); title("DENOISED IMAGE WITH AVERAGE FILTER - [k1 k1]");
subplot(4,2,8); imshow(Filtered_noise_image_3_4); title("DENOISED IMAGE WITH AVERAGE FILTER - [k2 k2]");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Exercise 4 - FLOATING PARTICLES                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
frame_4_1 = imread('initial_frame_1.jpg'); %Gathering Image (Resolution: 720x1080)
frame_4_2 = imread('initial_frame_2.jpg'); %Gathering Image (Resolution: 720x1080)
%k1, k2 remain the same
%Converting images to grayscale
gray_frame_4_1 = rgb2gray(frame_4_1);
gray_frame_4_2 = rgb2gray(frame_4_2);
%Applying noise reduction filter using median
Filtered_frame_4_1 = medfilt2(gray_frame_4_1, [k1 k1]);
Filtered_frame_4_2 = medfilt2(gray_frame_4_2, [k1 k1]);
%Applying h filter
h = [1 1 1; 1 -8 1; 1 1 1]; 
Filtered_frame_4_3 = imfilter(gray_frame_4_1, h);
Filtered_frame_4_4 = imfilter(gray_frame_4_2, h);
%Displaying and comparing images
figure(5); 
subplot(3,2,1); imshow(gray_frame_4_1); title("ORIGINAL FRAME 1");
subplot(3,2,2); imshow(gray_frame_4_2); title("ORIGINAL FRAME 2");
subplot(3,2,3); imshow(Filtered_frame_4_1); title("MEDIAN FILTERED FRAME 1");
subplot(3,2,4); imshow(Filtered_frame_4_2); title("MEDIAN FILTERED FRAME 2");
subplot(3,2,5); imshow(Filtered_frame_4_3); title("H FILTERED FRAME 1");
subplot(3,2,6); imshow(Filtered_frame_4_4); title("H FILTERED FRAME 2");