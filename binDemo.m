clc;
clear all;
close all;

%% Real Image
I = im2double(imread('4.bmp'));
if ndims(I) == 3
    I = rgb2gray(I);
end
I = I>0.5;
E = I;
p = 0.2;
rand('seed',sum(100*clock));
N = randi(1000,size(I,1),size(I,2));
N = N<=(p*1000);
E = xor(E,N);

%% Simulation Image
E = runReal(E);
m = max(E(:)); E = E./(m+(m==0));
figure;
subplot(1,2,1); imshow(I);
subplot(1,2,2); imshow(E);
