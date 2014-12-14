clc;
clear all;
close all;

%% Real Image
I = im2double(imread('images000.tiff'));
I = imresize(I,[129 129]); 
if ndims(I) == 3
    I = rgb2gray(I);
end

%% Simulation Image
I = im2double(imread('Images/CC.png'));
E = runReal(I);
m = max(E(:)); E = E./(m+(m==0));
figure;
subplot(1,2,1); imshow(I);
subplot(1,2,2); imshow(E);
