clc;
clear all;
close all;

I = im2double(imread('Images/CC.png'));
E = runReal(I);
m = max(E(:)); E = E./(m+(m==0));
figure;
subplot(1,2,1); imshow(I);
subplot(1,2,2); imshow(E);
