%% Parameters
function prm = getPrm()
    % Variables
    prm.minContrast = 9;
    prm.removeEpsilon = 0.248; 
    prm.nmsFact = 0.75;
    prm.maxTurn = 35;
    prm.maxNumOfEdges = 100;
    prm.addShift = true;
    
    % Constants
    prm.complexity = 15;
    prm.w = 2;
    prm.sigma = 0.1;
    prm.patchSize = 5;
    prm.fast = true;
    prm.run = true;
    prm.convertToSparse = false;
    prm.algoName = {'QP','TPT','OM','Canny','Lines','Sobel'};
    prm.algo = prm.algoName{1};
end