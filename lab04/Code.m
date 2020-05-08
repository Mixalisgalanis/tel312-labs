%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BASIC INFORMATION                                   %
% Course: Digital Image Processing - Lab 4            %
% Deadline: 07-05-2019                                %
% LAB31239720:  Pantelis Karamailis, 2016030040       %
%               Kostantinos Vlachos, 2016030042       %
%               Mixalis Galanis,     2016030036       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Clearing things up
clear all;
close all;
clearvars;

%Gathering Image and converting it to grayscale(Resolution: 510x343)
original_image = imread('image.jpg');
grayscale_image = rgb2gray(original_image); %Converting to grayscale
imwrite(grayscale_image, 'grayscale_image.jpg');
%Generating Gaussian Noise with 3 different noise
noises = [0.001 0.005 0.01];
noised_image_1 = imnoise(grayscale_image,'gaussian',0,noises(1));
noised_image_2 = imnoise(grayscale_image,'gaussian',0,noises(2));
noised_image_3 = imnoise(grayscale_image,'gaussian',0,noises(3));
%Blurring image with average filter
imshow(grayscale_image);
h = fspecial('average',[2 2]);
filtered_image_1 = imfilter(noised_image_1,h);
filtered_image_2 = imfilter(noised_image_2,h);
filtered_image_3 = imfilter(noised_image_3,h);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1 - WIENER RESTORATION FILTER                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculating Fourier Transforms of filter
[horizontal_pixels, vertical_pixels] = size(grayscale_image);
H = fft2(h,horizontal_pixels,vertical_pixels);
H_conj = conj(H);
%Calculating Fourier Transforms of Degraded Image
G_1 = fft2(filtered_image_1);
G_2 = fft2(filtered_image_2);
G_3 = fft2(filtered_image_3);
S_f = abs(fft2(grayscale_image)).^2;
%Calculating Fourier Transforms of Estimated Image
gamma = 10^7;
F_1_1 = (H_conj ./ (abs(H).^2 + gamma .* (noises(1) ./ S_f))).* G_1;  
F_1_2 = (H_conj ./ (abs(H).^2 + gamma .* (noises(2) ./ S_f))).* G_2; 
F_1_3 = (H_conj ./ (abs(H).^2 + gamma .* (noises(3) ./ S_f))).* G_3; 
%Estimated Image
IFT_1_1 = ifft2(F_1_1);
IFT_1_2 = ifft2(F_1_2);
IFT_1_3 = ifft2(F_1_3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2 - CONSTRAINED LEAST SQUARE RESTORATION FILTER     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P = [0 -1 0; -1 4 -1;0 -1 0];
M = horizontal_pixels + 2;
N = vertical_pixels + 2;
P_Extended = zeros(M,N);
P_Extended((M/2-1):(M/2+1),(N/2-1):(N/2+1)) = P;
P_F_T = fft2(P_Extended, M, N);
H = fft2(h,M,N);
H_conjugate = conj(H);

G_1 = fft2(filtered_image_1,M,N);
gamma_1 = 5;
norm_est_1 = norm((M-1)*(N-1)*(0 + noises(1)));
a_1 = norm_est_1 * 0.02;
while(1)
    F_1 = (H_conjugate./(abs(H).^2 + gamma_1*(abs(P_F_T).^2))).*G_1;
    f_1 = ifft2(F_1);
    R_1 = G_1 - H.*F_1;
    r_1 = ifft2(R_1);
    f_g_1 = norm(r_1);
    if(f_g_1<(norm_est_1 - a_1))
        gamma_1 = gamma_1 + 0.05 * gamma_1;
    end
    if (f_g_1>(norm_est_1 + a_1))
        gamma_1 = gamma_1 - 0.05 * gamma_1;
    end
    if((norm_est_1 - a_1) < f_g_1 && f_g_1 < (norm_est_1 + a_1))
        break;
    end
end

G_2 = fft2(filtered_image_2,M,N);
gamma_2 = 5;
norm_est_2 = norm((M-1)*(N-1)*(0 + noises(2)));
a_2 = norm_est_2 * 0.02;
while(1)
    F_2 = (H_conjugate./(abs(H).^2 + gamma_2*(abs(P_F_T).^2))).*G_2;
    f_2 = ifft2(F_2);
    R_2 = G_2 - H.*F_2;
    r_2 = ifft2(R_2);
    f_g_2 = norm(r_2);
    if(f_g_2<(norm_est_2 - a_2))
        gamma_2 = gamma_2 + 0.05 * gamma_2;
    end
    if (f_g_2>(norm_est_2 + a_2))
        gamma_2 = gamma_2 - 0.05 * gamma_2;
    end
    if((norm_est_2 - a_2) < f_g_2 && f_g_2 < (norm_est_2 + a_2))
        break;
    end
end

G_3 = fft2(filtered_image_3,M,N);
gamma_3 = 5;
norm_est_3 = norm((M-1)*(N-1)*(0 + noises(3)));
a_3 = norm_est_3 * 0.02;
while(1)
    F_3 = (H_conjugate./(abs(H).^2 + gamma_3*(abs(P_F_T).^2))).*G_3;
    f_3 = ifft2(F_3);
    R_3 = G_2 - H.*F_3;
    r_3 = ifft2(R_3);
    f_g_3 = norm(r_3);
    if(f_g_3<(norm_est_3 - a_3))
        gamma_3 = gamma_3 + 0.05 * gamma_3;
    end
    if (f_g_3>(norm_est_3 + a_3))
        gamma_3 = gamma_3 - 0.05 * gamma_3;
    end
    if((norm_est_3 - a_3) < f_g_3 && f_g_3 < (norm_est_3 + a_3))
        break;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUMMARY - COMPARING IMAGES AND CALCULATING MSE's    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Displaying and comparing noise images
subplot(4,1,1); imshow(grayscale_image); title("ORIGINAL GRAYSCALE IMAGE");
subplot(4,3,4); imshow(noised_image_1); title("AVERAGE FILTERED - GAUSSIAN NOISED IMAGE (0.001)");
subplot(4,3,5); imshow(noised_image_2); title("AVERAGE FILTERED - GAUSSIAN NOISED IMAGE (0.005)");
subplot(4,3,6); imshow(noised_image_3); title("AVERAGE FILTERED - GAUSSIAN NOISED IMAGE (0.01)");
subplot(4,3,7); imshow(uint8(IFT_1_1)); title("WIENER AVERAGE FILTERED - GAUSSIAN NOISED IMAGE (0.001)");
subplot(4,3,8); imshow(uint8(IFT_1_2)); title("WIENER AVERAGE FILTERED - GAUSSIAN NOISED IMAGE (0.005)");
subplot(4,3,9); imshow(uint8(IFT_1_3)); title("WIENER AVERAGE FILTERED - GAUSSIAN NOISED IMAGE (0.01)");
subplot(4,3,10); imshow(uint8(f_1)); title("CLSR AVERAGE FILTERED - GAUSSIAN NOISED IMAGE (0.001)");
subplot(4,3,11); imshow(uint8(f_2)); title("CLSR AVERAGE FILTERED - GAUSSIAN NOISED IMAGE (0.005)");
subplot(4,3,12); imshow(uint8(f_3)); title("CLSR AVERAGE FILTERED - GAUSSIAN NOISED IMAGE (0.01)");

% imwrite(noised_image_1, 'noised_image_0.001.jpg');
% imwrite(noised_image_2, 'noised_image_0.005.jpg');
% imwrite(noised_image_3, 'noised_image_0.01.jpg');
% imwrite(uint8(IFT_1_1), 'wiener_0.001.jpg');
% imwrite(uint8(IFT_1_2), 'wiener_0.005.jpg');
% imwrite(uint8(IFT_1_3), 'wiener_0.01.jpg');
% imwrite(uint8(f_1), 'clsr_0.001.jpg');
% imwrite(uint8(f_2), 'clsr_0.005.jpg');
% imwrite(uint8(f_3), 'clsr_0.01.jpg');



%Calculating Mean Squared Errors
[horizontal_pixels,vertical_pixels] = size(grayscale_image);
mse_1 = immse(grayscale_image, uint8(IFT_1_1));
mse_2 = immse(grayscale_image, uint8(IFT_1_2));
mse_3 = immse(grayscale_image, uint8(IFT_1_3));
mse_4 = immse(grayscale_image, uint8(f_1(2:end-1, 2:end-1)));
mse_5 = immse(grayscale_image, uint8(f_2(2:end-1, 2:end-1)));
mse_6 = immse(grayscale_image, uint8(f_3(2:end-1, 2:end-1)));