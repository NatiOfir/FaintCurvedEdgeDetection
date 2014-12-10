function [E,I] = run( I,prm,disp )
    clc;
    close all;
    
    if ~exist('prm','var')
        prm = getPrm();
    end
    
    if ~exist('disp','var')
        disp = false;
    end
    
    I = im2double(I);
    res = runIm(I,prm);
    E = res.E;
    
    if prm.addShift
        Inew = zeros(size(I));
        Is = I(3:end,3:end);
        Inew(1:end-2,1:end-2) = Is;
        Inew(end-1:end,1:end) = I(1:2,1:end);
        Inew(1:end,end-1:end) = I(1:end,1:2);
        res = runIm(Inew,prm);
        Enew = res.E;
        E(3:end,3:end) = max(E(3:end,3:end),Enew(1:end-2,1:end-2));
    end
    
    if disp
        E = E./max(E(:));
        [a,name,c] = fileparts(imgStr);
        imwrite(E,sprintf('%s_res.png',name));
        figure;
        subplot(1,2,1);imshow(I);
        subplot(1,2,2);imshow(E);
    end
end

