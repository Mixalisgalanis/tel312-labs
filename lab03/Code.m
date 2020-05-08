%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BASIC INFORMATION                                   %
% Course: Digital Image Processing - Lab 3            %
% Deadline: 09-04-2019                                %
% LAB31239720:  Pantelis Karamailis, 2016030040       %
%               Kostantinos Vlachos, 2016030042       %
%               Mixalis Galanis,     2016030036       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Clearing things up
close all;
clearvars;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Exercise 1 - FOURIER TRANSFORMATION OF TOOLS IMAGE   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
original_image_1_1 = imread('tools.bmp'); %Gathering Image (Resolution: 510x343)
%Applying Fourier Transformation to the image without shifting
F_1_1 = abs(fft2(original_image_1_1));
%Applying Fourier Transformation to the image with shifting
F_1_2 = abs(fftshift(fft2(original_image_1_1)));
%Applying Log Filters
c_1_2 = 0.04; % 0.02 < c < 0.07 works best for most spectrum detail in grayscale
Fourier_image_1_1 = c_1_2 * log(1 + F_1_1);
Fourier_image_1_2 = c_1_2 * log(1 + F_1_2);
%Using Color Map to Display Fourier Transformation
colormap(jet());
c_1_3 = 4; % 3 < c < 4.5 works best for most spectrum detail with color map
Fourier_image_1_3 = c_1_3 * log(1 + F_1_2);
%Displaying and comparing images
figure(1);
subplot(2,2,1); imshow(original_image_1_1); title("ORIGINAL IMAGE");
subplot(2,2,2); imshow(Fourier_image_1_1); title("FOURIER TRANSFORM OF IMAGE (NO SHIFT)");
subplot(2,2,3); imshow(Fourier_image_1_2); title("FOURIER TRANSFORM OF IMAGE (SHIFT, GRAY)");
subplot(2,2,4); image(Fourier_image_1_3); title("FOURIER TRANSFORM OF IMAGE (SHIFT, RGB)");
%imwrite(Fourier_image_1_1, 'Fourier_image_1_1.jpg');
%imwrite(Fourier_image_1_2, 'Fourier_image_1_2.jpg');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Exercise 2 - FOURIER TRANSFORMATION OF IMAGE PATTERS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
original_image_2_1 = (imread('horizontal_lines_pattern.jpg'));
original_image_2_2 = imread('vertical_lines_pattern.jpg');
%Applying Fourier Transformation to the image with shifting
F_2_1 = abs(fftshift(fft2(original_image_2_1)));
F_2_2 = abs(fftshift(fft2(original_image_2_2)));
%Applying Log Filters
c_2_1 = 0.04; % 0.02 < c < 0.07 works best for most spectrum detail in grayscale
c_2_2 = 0.02; % 0.02 < c < 0.07 works best for most spectrum detail in grayscale
Fourier_image_2_1 = c_2_1 * log(1 + F_2_1);
Fourier_image_2_2 = c_2_2 * log(1 + F_2_2);
%Displaying and comparing images
figure(2);
subplot(1,2,1); imshow(original_image_2_1); title("HORIZONTAL LINES (ORIGINAL)");
subplot(1,2,2); imshow(Fourier_image_2_1); title("HORIZONTAL LINES (FOURIER TRANSFORM)");
figure(3);
subplot(1,2,1); imshow(original_image_2_2); title("VERTICAL LINES (ORIGINAL)");
subplot(1,2,2); imshow(Fourier_image_2_2); title("VERTICAL LINES (FOURIER TRANSFORM)");
%imwrite(Fourier_image_2_1, 'Fourier_image_2_1.jpg');
%imwrite(Fourier_image_2_2, 'Fourier_image_2_2.jpg');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Exercise 3 - FOURIER TRANSFORMATION FOR BLACK IMAGE  %
%              WITH WHITE BOX                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Creating image size
image_resolution = 256;
diameters = [10 30];
%Analyzing Pixel Data
for m = 1 : length(diameters) %For every diameter
    %Creating image and target area
    image = zeros(image_resolution);
    target_space_start = (image_resolution - diameters(m))/2;
    target_space_end = (image_resolution + diameters(m))/2;
    %Filling White Box
    for i = 1 : length(image) %For every pixel in horizontal axis
        for j = 1 : length(image) %For every pixel in vertical axis
            if (i >= target_space_start && i <= target_space_end && j >= target_space_start && j <= target_space_end)
                image(i,j) = 255;
            end
        end 
    end
    %Applying Fourier Transformation & Log Filters
    F_3_1 = abs(fftshift(fft2(image)));
    c_3_1 = 0.04; % 0.02 < c < 0.07 works best for most spectrum detail in grayscale
    Fourier_image_3_1 =  c_3_1 * log(1 + F_3_1);
    %Displaying and comparing images
    figure(4);
    subplot(2,2,2*m - 1); imshow(image); title([num2str(diameters(m)) 'x' num2str(diameters(m)) ' diameter box (ORIGINAL)']);
    subplot(2,2,2*m); imshow(Fourier_image_3_1); title([num2str(diameters(m)) 'x' num2str(diameters(m)) ' diameter box (FOURIER TRANSFORM)']);
    %imwrite(image, ['Fourier_image_3_1_' num2str(2*m - 1) '.jpg']);
    %imwrite(Fourier_image_3_1, ['Fourier_image_3_1_' num2str(2*m) '.jpg']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Exercise 4 - FOURIER TRANSFORMATION FOR BLACK IMAGE  %
%              WITH ROTATED WHITE BOX                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Creating image size
image_resolution = 256;
diameters = [10 30];
%Analyzing Pixel Data
for m = 1 : length(diameters) %For every diameter
    %Creating image and target area
    image = zeros(image_resolution);
    target_space_start = (image_resolution - diameters(m))/2;
    target_space_end = (image_resolution + diameters(m))/2;
    %Filling White Box
    for i = 1 : length(image) %For every pixel in horizontal axis
        for j = 1 : length(image) %For every pixel in vertical axis
            if (i >= target_space_start && i <= target_space_end && j >= target_space_start && j <= target_space_end)
                image(i,j) = 255;
            end
        end 
    end
    %Rotating image first
    image = imrotate(image, 45);
    %Applying Fourier Transformation & Log Filters
    F_4_1 = abs(fftshift(fft2(image)));
    c_4_1 = 0.04; % 0.02 < c < 0.07 works best for most spectrum detail in grayscale
    Fourier_image_4_1 =  c_4_1 * log(1 + F_4_1);
    %Displaying and comparing images
    figure(5);
    subplot(2,2,2*m - 1); imshow(image); title([num2str(diameters(m)) 'x' num2str(diameters(m)) ' diameter box rotated (ORIGINAL)']);
    subplot(2,2,2*m); imshow(Fourier_image_4_1); title([num2str(diameters(m)) 'x' num2str(diameters(m)) ' diameter box rotated (FOURIER TRANSFORM)']);
    %imwrite(image, ['Fourier_image_4_1_' num2str(2*m - 1) '.jpg']);
    %imwrite(Fourier_image_4_1, ['Fourier_image_4_1_' num2str(2*m) '.jpg']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Exercise 5 - FOURIER TRANSFORMATION FOR BLACK IMAGE  %
%              WITH RELOCATED WHITE BOX                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Creating image size
image_resolution = 256;
diameters = [10 30];
%Generating random offsets
random_offsets = zeros(1,2);
for k = 1 : 2 % 1 random offset for X axis, 1 random offset for Y axis
    random_offsets(k) = (rand - 0.5) * (image_resolution - max(diameters));
end
%Analyzing Pixel Data
for m = 1 : length(diameters) %For every diameter
    %Creating image and target area with random offset
    image = zeros(image_resolution);
    target_space_start_X = (image_resolution - diameters(m))/2 + random_offsets(1);
    target_space_end_X = (image_resolution + diameters(m))/2 + random_offsets(1);
    target_space_start_Y = (image_resolution - diameters(m))/2 + random_offsets(2);
    target_space_end_Y = (image_resolution + diameters(m))/2 + random_offsets(2);
    %Filling White Box
    for i = 1 : length(image) %For every pixel in horizontal axis
        for j = 1 : length(image) %For every pixel in vertical axis
            if (i >= target_space_start_X && i <= target_space_end_X && j >= target_space_start_Y && j <= target_space_end_Y)
                image(i,j) = 255;
            end
        end 
    end
    %Applying Fourier Transformation & Log Filters
    F_5_1 = abs(fftshift(fft2(image)));
    c_5_1 = 0.04; % 0.02 < c < 0.07 works best for most spectrum detail in grayscale
    Fourier_image_5_1 =  c_5_1 * log(1 + F_3_1);
    %Displaying and comparing images
    figure(6);
    subplot(2,2,2*m - 1); imshow(image); title([num2str(diameters(m)) 'x' num2str(diameters(m)) ' diameter box rotated (ORIGINAL)']);
    subplot(2,2,2*m); imshow(Fourier_image_5_1); title([num2str(diameters(m)) 'x' num2str(diameters(m)) ' diameter box rotated (FOURIER TRANSFORM)']);
    %imwrite(image, ['Fourier_image_5_1_' num2str(2*m - 1) '.jpg']);
    %imwrite(Fourier_image_5_1, ['Fourier_image_5_1_' num2str(2*m) '.jpg']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Exercise 6 - DIFFERENCE BETWEEN ORIGINAL & EQUALIZED %
%              IMAGE                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
original_image_6_1 = imread('tools.bmp'); %Gathering Image (Resolution: 510x343)
equalized_image_6_2 = histeq(original_image_6_1);
%Applying Fourier Transformations to the images with shifting
F_6_1 = abs(fftshift(fft2(original_image_6_1)));
F_6_2 = abs(fftshift(fft2(equalized_image_6_2)));
%Applying Log Filters
c_6_1 = 0.04; % 0.02 < c < 0.07 works best for most spectrum detail in grayscale
c_6_2 = 0.04; % 0.02 < c < 0.07 works best for most spectrum detail in grayscale
Fourier_image_6_1 = c_6_1 * log(1 + F_6_1);
Fourier_image_6_2 = c_6_2 * log(1 + F_6_2);
%Displaying and comparing images
figure(7);
subplot(2,2,1); imshow(original_image_6_1); title("ORIGINAL IMAGE");
subplot(2,2,2); imshow(equalized_image_6_2); title("EQUALIZED IMAGE");
subplot(2,2,3); imshow(Fourier_image_6_1); title("FOURIER TRANSFORM OF ORIGINAL IMAGE");
subplot(2,2,4); imshow(Fourier_image_6_2); title("FOURIER TRANSFORM OF EQUALIZED IMAGE");
%imwrite(equalized_image_6_2, ['equalized_image_6_1.jpg']);
%imwrite(Fourier_image_6_1, ['Fourier_image_6_2.jpg']);
%imwrite(Fourier_image_6_2, ['Fourier_image_6_3.jpg']);