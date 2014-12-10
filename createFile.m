clc;
clear;
close all;

prm = getPrm();
prm.run = false;
prm.fast = false;

if prm.complexity == 2
    d1 = 'Files2'; d2 = 'Files2_2';
elseif prm.complexity == 15
    d1 = 'Files15'; d2 = 'Files15_2';
end
if exist(d1,'dir')
    %rmdir(d1,'s');
end
if exist(d2,'dir')
    %rmdir(d2,'s');
end
profile on;
J = 10;
J0 = 9;
for j=J0:J
    I = zeros(2^j+1);
    runIm(I,prm);
    clear functions;
end

prm.fast = true;

for j=J0:J
    I = zeros(2^j+1);
    runIm(I,prm);
    clear funcions;
end
profile viewer;
profile off;