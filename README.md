FaintCurvedEdgeDetection
========================

Detector of faint edges in noisy images for Matlab.
This detector is relatively slow, I suggest to use one of my two other solutions for faint edge detection.

Usage in Matlab:

>> I  = im2double(imread('img.png'));
>> I = runIm(I);

See demo.m for example.

Important Files:
getPrm.m - Params of the algorithm
