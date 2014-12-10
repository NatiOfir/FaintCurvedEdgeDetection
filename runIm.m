%% Main Function

function resOut = runIm(I,prm)
    disp('Build Quad Pyramid');
    tic;
    global res;

    if ~exist('prm','var')
        prm = getPrm();
    end
    res.prm = prm;
    [n,m] = size(I);
    j = log2(n-1);
    
    if n ~= m
        disp('Non square image');
    elseif j-floor(j) ~=0
        disp('Image size is not a power of 2 plus 1');
    else 
        bot = extractBottomLevelData(5,prm);
        
        g = meshgrid(1:n);
        S = sub2ind([n,n],g',g);
        res.S = S;
        %disp(S);
        N = n^2;
        res.N = N;
        sparseFlag = false;
        fClass = 'double';
        iClass = 'uint16';
        
        if prm.run
            sparseFlag = true;
            res.R = nan(N,N,fClass);
            res.L = nan(N,N,fClass);
            res.C = nan(N,N,fClass);
            res.Ang0 = nan(N,N,fClass);
            res.Ang1 = nan(N,N,fClass);
            res.minC = nan(N,N,fClass);
            res.maxC = nan(N,N,fClass);
            res.S0 = zeros(N,N,iClass);
            if prm.complexity == 2
                res.S1 = zeros(N,N,iClass);
            end
            res.pixels = zeros(N,N,prm.patchSize,iClass);
        else
            res.R = sparse(N,N);
            res.L = sparse(N,N);
            res.C = sparse(N,N);
            res.Ang0 = sparse(N,N);
            res.Ang1 = sparse(N,N);
            res.minC = sparse(N,N);
            res.maxC = sparse(N,N);
            res.S0 = sparse(N,N);
            if prm.complexity == 2
                res.S1 = sparse(N,N);
            end
        end
        
        if prm.fast
            QPFast(I,S,bot,N,prm);
        else
            QP(I,S,bot,N,prm);
        end
        
        if prm.minContrast>0 && prm.run
            minTest = ((res.C>0 & res.minC>=(res.C/2)) | (res.C<0 & res.maxC<=(res.C/2)));
            res.C(~minTest) = nan;
        end
        toc;
        
        disp('Create Edge Image');
        tic;
        if prm.run
            clear functions;
            getScores(prm);
        end
        toc;
        
        disp('Convert to Sparse');
        tic;
        if sparseFlag && prm.convertToSparse
            res.R = convertToSparse(res.R);
            res.L = convertToSparse(res.L);
            res.C = convertToSparse(res.C);
            res.Ang0 = convertToSparse(res.Ang0);
            res.Ang1 = convertToSparse(res.Ang1);
            res.minC = convertToSparse(res.minC);
            res.maxC = convertToSparse(res.maxC);
            res.S0 = convertToSparse(res.S0);
            if prm.complexity == 2
                res.S1 = convertToSparse(res.S1);
            end
        end
        toc;
    end
    
    resOut = res;
    clear res;
end

function D = convertToSparse(M)
    M(isnan(M)) = 0;
    D = sparse(double(M));
end

%% Quad Pyramid Recursive

function QPFast(I,S,bot,N,prm)
    global res;
    
    n = sqrt(N);
    w = prm.w;
    
    if prm.run
        step = prm.patchSize-1;
        Ipad = padarray(I,[w w],'symmetric','both');
        Ipatches = im2col(Ipad,[2*w+step+1 2*w+step+1],'sliding');
        Spatches = im2col(S,[step+1 step+1],'sliding');
        
        xInd = floor(step/2)+1:step:n-floor(step/2);
        xInd = xInd-floor(step/2);
        [X,Y] = meshgrid(xInd,xInd);
        Ind = sub2ind([n-step n-step],X,Y);
        Ind = Ind(:)';
        
        Ipatches = Ipatches(:,Ind);
        Spatches = Spatches(:,Ind);
        
        getBottomLevelFast(Ipatches,Spatches,bot,prm);
    end
    
    maxJ = log2(n-1)-1;
    for j = 2:maxJ
        m = 2^j+1;
        prefix = sprintf('merge%d_%d_%d',prm.complexity,n,m);
        fileName = sprintf('Files%d_2/%s.mat',prm.complexity,prefix);
                
        if ~exist(fileName,'file')
            dirName = sprintf('Files%d_2',prm.complexity);
            if ~exist(dirName,'dir')
                mkdir(dirName);
            end
            createMergeFile(prefix,fileName,prm);
        end
        if prm.run
            load(fileName)
            disp(fileName)
            mergeSquaresFast(tableSingle,prm);
            clear functions;
            clear tableSingle;
            mergeSquaresFast(tableDouble,prm);
            clear function;
            clear tableDouble;
        end
    end
end

function QP(I,S,bot,N,prm)    
    [n,m] = size(I);
    
    if n == 5 
        if prm.run
            getBottomLevel(I,S,bot,prm);
        end
    else
        mid = ceil(n/2);

        c = cell(2,2);
        [c{1,1}.I,c{1,1}.S] = subIm(I,S,1,1,mid,mid);
        [c{1,2}.I,c{1,2}.S] = subIm(I,S,mid,1,n,mid);
        [c{2,1}.I,c{2,1}.S] = subIm(I,S,1,mid,mid,n);
        [c{2,2}.I,c{2,2}.S] = subIm(I,S,mid,mid,n,n);
        
        for i=1:2
            for j=1:2
                QP(c{i,j}.I,c{i,j}.S,bot,N,prm);
            end
        end
        
        if prm.complexity == 2
            mergeSquares2(c{1,1}.S,c{2,1}.S,c{1,2}.S,c{2,2}.S,N,prm);
        elseif prm.complexity == 15
            mergeSquares15(c{1,1}.S,c{2,1}.S,c{1,2}.S,c{2,2}.S,N,prm);
        end
    end
end

%% Pyramid Building

function mergeSquares2(S11,S21,S12,S22,N,prm)

    global res;
    
    [n n] = size(S11);
      
    dirName = sprintf('Files%d',prm.complexity);
    
    if ~exist(dirName,'dir');
        mkdir(dirName);
    end
    
    fileName = sprintf('Files%d/merge%d_%d_%d_%d.mat',prm.complexity,prm.complexity,sqrt(N),n,S11(1,1));
    
    if exist(fileName,'file')
        if prm.run
            load(fileName)
            mergeSquaresFast(tableSingle,prm);
            mergeSquaresFast(tableDouble,prm);
            clear function;
            clear tableSingle;
            clear tableDouble;
        end
        return;
    end
    
    disp(fileName);
    disp('Do not exist');
    
    [tableSingleFinal,tableDoubleFinal1,tableDoubleFinal2] = merge2(S11,S21,S12,S22); 
    tableDoubleFinal = [tableDoubleFinal1;tableDoubleFinal2];

    clear functions;
        
    rSize = size(res.L);
    
    ind0 = tableSingleFinal(:,1);
    ind1 = tableSingleFinal(:,2);
    s0 = tableSingleFinal(:,3);
    %angle = getAngle(ind0,ind1,s0);
    ind0s0 = sub2ind(rSize,ind0,s0);
    s0ind1 = sub2ind(rSize,s0,ind1);
    ind01 = sub2ind(rSize,ind0,ind1);
    ind10 = sub2ind(rSize,ind1,ind0);
    tableSingleFinal = [s0,ind0s0,s0ind1,ind01,ind10];
 
    ind0 = tableDoubleFinal(:,1);
    ind1 = tableDoubleFinal(:,2);
    s0 = tableDoubleFinal(:,3);
    %angle = getAngle(ind0,ind1,s0);
    s1 = tableDoubleFinal(:,4);
    ind0s0 = sub2ind(rSize,ind0,s0);
    s0s1 = sub2ind(rSize,s0,s1);
    s1ind1 = sub2ind(rSize,s1,ind1);
    ind01 = sub2ind(rSize,ind0,ind1);
    ind10 = sub2ind(rSize,ind1,ind0);
    tableDoubleFinal = [s0,s1,ind0s0,s0s1,s1ind1,ind01,ind10];
    
    s.tableSingle = tableSingleFinal;
    s.tableDouble = tableDoubleFinal;
        
    save(fileName,'-struct','s','-v7.3');    
end

function mergeSquares15(S11,S21,S12,S22,N,prm)

    global res;
    
    [n n] = size(S11);
    
    dirName = sprintf('Files%d',prm.complexity);
    
    if ~exist(dirName,'dir');
        mkdir(dirName);
    end
    
    fileName = sprintf('Files%d/merge%d_%d_%d_%d.mat',prm.complexity,prm.complexity,sqrt(N),n,S11(1,1));
        
    if exist(fileName,'file')
        if prm.run
            load(fileName)
            mergeSquaresFast(tableSingle,prm);
            mergeSquaresFast(tableDouble,prm);
            clear function;
            clear tableSingle;
            clear tableDouble;
        end
        return;
    end
    
    disp(fileName);
    disp('Do not exist');
    
    [tableSingleFinal,tableDoubleFinal] = merge15(S11,S21,S12,S22);
    clear functions;
    
    rSize = size(res.L);
    
    ind0 = tableSingleFinal(:,1);
    ind1 = tableSingleFinal(:,2);
    s0 = tableSingleFinal(:,3);
    %angle = getAngle(ind0,ind1,s0);
    
    ind0s0 = sub2ind(rSize,ind0,s0);
    s0ind1 = sub2ind(rSize,s0,ind1);
    ind01 = sub2ind(rSize,ind0,ind1);
    ind10 = sub2ind(rSize,ind1,ind0);
    tableSingleFinal = [s0,ind0s0,s0ind1,ind01,ind10];
 
    ind0 = tableDoubleFinal(:,1);
    ind1 = tableDoubleFinal(:,2);
    s0 = tableDoubleFinal(:,3);
    %angle = getAngle(ind0,ind1,s0);
    ind0s0 = sub2ind(rSize,ind0,s0);
    s0ind1 = sub2ind(rSize,s0,ind1);
    ind01 = sub2ind(rSize,ind0,ind1);
    ind10 = sub2ind(rSize,ind1,ind0);
    tableDoubleFinal = [s0,ind0s0,s0ind1,ind01,ind10];
    
    s.tableSingle = tableSingleFinal;
    s.tableDouble = tableDoubleFinal;
        
    save(fileName,'-struct','s','-v7.3');
end

function mergeSquaresFast(table,prm)
    global res;
    
    single = (size(table,2) == 5);
    
    if single
        s0 = table(:,1);
        ind0s0 = table(:,2);
        s0ind1 = table(:,3);
        ind01 = table(:,4);
        ind10 = table(:,5);
             
        validLen = res.L(ind0s0)>=0 & res.L(s0ind1)>=0;
        len = res.L(ind0s0)+res.L(s0ind1);
        resp = res.R(ind0s0)+res.R(s0ind1);
        stitchAng = mod(res.Ang1(ind0s0)-res.Ang0(s0ind1),360);
        validAng = stitchAng<=prm.maxTurn | 360-stitchAng<=prm.maxTurn;
        minC = min(res.minC(ind0s0),res.minC(s0ind1));
        maxC = max(res.maxC(ind0s0),res.maxC(s0ind1));
        ang0 = res.Ang0(ind0s0);
        ang1 = res.Ang1(s0ind1);
    else
        s0 = table(:,1);
        s1 = table(:,2);
        ind0s0 = table(:,3);
        s0s1 = table(:,4);
        s1ind1 = table(:,5);
        ind01 = table(:,6);
        ind10 = table(:,7);

        validLen = res.L(ind0s0)>=0 & res.L(s0s1)>=0 & res.L(s1ind1)>=0;
        len = res.L(ind0s0)+res.L(s0s1)+res.L(s1ind1);
        resp = res.R(ind0s0)+res.R(s0s1)+res.R(s1ind1);
        stitchAng1 = mod(res.Ang1(ind0s0)-res.Ang0(s0s1),360);
        validAng1 = stitchAng1<=prm.maxTurn | 360-stitchAng1<=prm.maxTurn;
        stitchAng2 = mod(res.Ang1(s0s1)-res.Ang0(s1ind1),360);
        validAng2 = stitchAng2<=prm.maxTurn | 360-stitchAng2<=prm.maxTurn;
        validAng = validAng1 & validAng2;
        minC = min(min(res.minC(ind0s0),res.minC(s0s1)),res.minC(s1ind1));
        maxC = max(max(res.maxC(ind0s0),res.maxC(s0s1)),res.maxC(s1ind1));
        ang0 = res.Ang0(ind0s0);
        ang1 = res.Ang1(s1ind1);
    end
    
    if prm.w == 0
        con = resp./len;
    else
        con = resp./len/prm.w/2;
    end
    con(~validLen) = nan;
    resp(~validLen) = nan;
    
    minLen = len<=prm.minContrast;
    minC(minLen) = con(minLen);
    maxC(minLen) = con(minLen);
        
    newRes = (abs(con) > abs(res.C(ind01))) | isnan(res.C(ind01));
    valid = (validLen>0 & newRes & validAng);
    
    if single
        data = [con len resp minC maxC ind01 ind10 s0 ang0 ang1];
    else
        data = [con len resp minC maxC ind01 ind10 s0 ang0 ang1 s1];
    end
    data = data(valid,:);
    scores = abs(data(:,1))-threshold(data(:,2));
    [Y,I] = sort(scores,'ascend');
    data = data(I,:);

    con = data(:,1);
    len = data(:,2);
    resp = data(:,3);
    minC = data(:,4);
    maxC = data(:,5);
    ind01 = data(:,6);
    ind10 = data(:,7);
    s0 = data(:,8);
    ang0 = data(:,9);
    ang1 = data(:,10);
    if ~single
        s1 = data(:,11);
    end
    
    res.R(ind01) = resp;
    res.R(ind10) = -resp;
    res.L(ind01) = len;
    res.L(ind10) = len;
    res.C(ind01) = con;
    res.C(ind10) = -con;
    res.Ang0(ind01) = ang0;
    res.Ang1(ind01) = ang1;
    res.Ang0(ind10) = mod(ang1+180,360);
    res.Ang1(ind10) = mod(ang0+180,360);
    res.minC(ind01) = minC;
    res.minC(ind10) = -maxC;
    res.maxC(ind01) = maxC;
    res.maxC(ind10) = -minC;
    
    if single
        res.S0(ind01) = s0;
        res.S0(ind10) = s0;
    else
        res.S0(ind01) = s0;
        res.S1(ind01) = s1;
        res.S0(ind10) = s1;
        res.S1(ind10) = s0;
    end
end

function getBottomLevelFast(Ipatches,Spatches,bot,prm)

    global res;
    
    if prm.w == 0
        resp = (bot.lineVec)'*Ipatches;
    else
        resp = (bot.leftVec-bot.rightVec)'*Ipatches;
    end
        
    ind0 = Spatches(bot.p0,:);
    ind1 = Spatches(bot.p1,:);
    angle = getAngle(ind0,ind1,sqrt(res.N));
    indices = bot.indices(:);
    isZero = (indices == 0);
    indices(isZero) = 1;
    lineIndices = Spatches(indices,:);
    lineIndices(isZero,:) = 0;
    
    len = repmat(bot.lengthVec',[1 size(Spatches,2) ]);
    
    rSize = size(res.L);
    ind01 = sub2ind(rSize,ind0,ind1);
    ind10 = sub2ind(rSize,ind1,ind0);
    
    N = size(res.L,1);

    for cord = 1:prm.patchSize
        curCord = lineIndices(cord:(prm.patchSize):end,:);
        temp = zeros(N,N);
        temp(ind01) = curCord;
        temp(ind10) = curCord;
        res.pixels(:,:,cord) = temp;
    end

    if prm.w == 0
        con = resp./len;
    else
        con = resp./len/prm.w/2;
    end
    
    con(len == 0 | isnan(len)) = nan;
    
    if prm.removeEpsilon>0
        bad = abs(con)<(prm.removeEpsilon*prm.sigma);
        con(bad) = nan;
        len(bad) = nan;
        resp(bad) = nan;
    end
    
    res.R(ind01) = resp;
    res.R(ind10) = -resp;
    res.L(ind01) = len;
    res.L(ind10) = len;
    res.C(ind01) = con;
    res.C(ind10) = -con;
    res.Ang0(ind01) = angle;
    res.Ang1(ind01) = angle;
    angle2 = mod(angle+180,360);
    res.Ang0(ind10) = angle2;
    res.Ang1(ind10) = angle2;
end

function getBottomLevel(I,S,bot,prm)

    global res;

    [n n] = size(I);
    I = reshape(I,[n^2 1]);
    
    if prm.w == 0
        resp = bsxfun(@times,bot.lineVec,I);
    else
        resp = bsxfun(@times,(bot.leftVec-bot.rightVec),I);
    end
    
    resp = sum(resp);
    
    ind0 = S(bot.p0);
    ind1 = S(bot.p1);
    len = bot.lengthVec;
    
    rSize = size(res.L);
    ind01 = sub2ind(rSize,ind0,ind1);
    ind10 = sub2ind(rSize,ind1,ind0);
    
    res.R(ind01) = resp;
    res.R(ind10) = -resp;
    
    if prm.w == 0
        con = resp./len;
    else
        con = resp./len/prm.w/2;
    end
    
    con(len == 0 | isnan(len)) = nan;
    
    if prm.removeEpsilon>0
        bad = abs(con)<(prm.removeEpsilon*prm.sigma);
        con(bad) = nan;
        len(bad) = nan;
        resp(bad) = nan;
    end
    
    res.L(ind01) = len;
    res.L(ind10) = len;
    res.C(ind01) = con;
    res.C(ind10) = -con;
end

%% Preprocessing for the Bottom Level Filters

function bot = extractBottomLevelData(n,prm)

    w = prm.w;
    nBig = n+2*w;
    
    fileName = sprintf('Files/bot%d_%d.mat',n,w);
    
    if exist(fileName,'file')
        load(fileName);
        return;
    end
    
    pairs = nchoosek(4,2)*n^2;
    lineVec = zeros(n^2,pairs);
    leftVec = zeros(nBig^2,pairs);
    rightVec = zeros(nBig^2,pairs);
    lengthVec = zeros(1,pairs);
    p0 = zeros(1,pairs);
    p1 = zeros(1,pairs);
    indices = zeros(prm.patchSize,pairs);
    
    index = 1;
    
    for e0 = 1:3
        for e1 = e0+1:4
            for v0 = 1:n
                [x0,y0] = getVerticesFromPatchIndices(e0,v0,n);
                for v1 = 1:n
                    [x1,y1] = getVerticesFromPatchIndices(e1,v1,n);
                    p0(index) = sub2ind([n,n],x0,y0);
                    p1(index) = sub2ind([n,n],x1,y1);
                    [P,L] = getLine(n,x0,y0,x1,y1);
                    curIndices = find(P);
                    indices(1:length(curIndices),index) = curIndices;
                    
                    line = P;
                    left = zeros(nBig);
                    right = zeros(nBig);
                    
                    dx = x1-x0;
                    dy = y1-y0;
                    
                    if abs(dx)>abs(dy)
                        dy = -sign(dx);
                        dx = 0;
                    else
                        dx = sign(dy);
                        dy = 0;
                    end
                    dx = sign(dx);
                    dy = sign(dy);
                    
                    for k=1:w
                        x = w+dx*k+1; y = w+dy*k+1;
                        right(x:x+n-1,y:y+n-1) = right(x:x+n-1,y:y+n-1)+line;
                        x = w-dx*k+1; y = w-dy*k+1;
                        left(x:x+n-1,y:y+n-1) = left(x:x+n-1,y:y+n-1)+line;
                    end
                    
                    right(right>1) = 1;
                    left(left>1) = 1;
                    
                    %figure, imshow(imresize((right-left+1)/2,30));
                    %disp(L);
                    %disp(right-left);
                    %fprintf('\n(%d,%d)\n',x0,y0);
                    
                    lineVec(:,index) = reshape(line,[n^2 1]);
                    leftVec(:,index) = reshape(left,[nBig^2 1]);
                    rightVec(:,index) = reshape(right,[nBig^2 1]);
                    lengthVec(index) = L;
                    index = index+1;
                end
            end
        end
    end
    
    bot.lineVec = lineVec;
    bot.leftVec = leftVec;
    bot.rightVec = rightVec;
    bot.lengthVec = lengthVec;
    bot.pairs = pairs;
    bot.p0 = p0;
    bot.p1 = p1;
    bot.indices = indices;
    s.bot = bot;
    save(fileName,'-struct','s','-v7.3');
end

function [P,L] = getLine(n,x0,y0,x1,y1)
    P = zeros(n);
    dx = abs(x1-x0);
    dy = abs(y1-y0);
    L = max(dx,dy);

    if x0 < x1 
        sx = 1;
    else
        sx = -1;
    end

    if y0 < y1
        sy = 1;
    else
        sy = -1;
    end

    err = dx-dy;
    first = true;
    
    while 1>0
        
        P(x0,y0) = 1;
        if first 
            P(x0,y0) = 0.5;
            first = false;
        end
            
        if x0 == x1 && y0 == y1
            P(x0,y0) = 0.5;
            break;
        end

        e2 = 2*err;

        if e2 > -dy 
            err = err - dy;
            x0 = x0 + sx;
        end
        if e2 < dx
            err = err + dx;
            y0 = y0 + sy;
        end
    end
end

function [x,y] = getVerticesFromPatchIndices(e,v,n)
    switch e
        case 1
            x = v;
            y = 1;
        case 2
            x = 1;
            y = v;
        case 3
            x = v;
            y = n;
        case 4
            x = n;
            y = v;
        otherwise
            x = -1;
            y = -1;
    end
end

%% Post Processing

function getScores(prm)
    
    global res;
    
    N = size(res.L,1);
    n = sqrt(N);
    res.E = zeros(n);
    selected = false(n);
    edgeCounter = 0;
    
    L = res.L;
    T = threshold(L);
    SC = abs(res.C)-T;
    LC = res.L;
    IN = ones(size(SC));
    IN = tril(IN);
    SC = tril(SC);
    LC = tril(LC);
    SC(IN == 0) = nan;
    LC(IN == 0) = nan;
    
    Ind = 1:(N^2);
    SC = SC(:); Ind = Ind(:); LC = LC(:);
    edge = SC>0 & ~isnan(SC);
    SC = SC(edge);
    LC = LC(edge);
    Ind = Ind(edge);
    
    [SC,I] = sort(SC,'descend');
    Ind = Ind(I);
    %LC = LC(I);
    %[LC,I] = sort(LC,'descend');
    %Ind = Ind(I);
    %SC = SC(I);
    counter = 0;
    for i = 1:length(SC)
        score = SC(i);
        curInd = Ind(i);
        tSize = [N,N];
        [ind0,ind1] = ind2sub(tSize,curInd);
        E = false(n);
        E = addEdge(ind0,ind1,tSize,E,prm,1);

        if E == -1
            continue;
        end
        Ebin = E>0;
        %fprintf('L = %d',sum(Ebin(:)));
        if false;
            M = zeros(size(E));
            M(ind0) = 1;
            M(ind1) = 1;
            subplot(1,2,1);imshow(E);
            subplot(1,2,2);imshow(M);
        end
        
        if prm.nmsFact == 0
            res.E = max(res.E,E*score);
        else
            curI = E;
            imSize = size(curI);
            curIdialate = curI | [zeros(1,imSize(2)) ; curI(1:end-1,1:end)]|[curI(2:end,1:end);zeros(1,imSize(2))]|[curI(1:end,2:end) zeros(imSize(1),1)]|[zeros(imSize(1),1) curI(1:end,1:end-1)];
            L = sum(curI(:));
            coor = curIdialate & selected;
            nmsScore = sum(coor(:))/L;
            if nmsScore < prm.nmsFact
                counter = counter+1;
                %figure,imshow(curI);
                edgeCounter = edgeCounter+1;
                selected = selected | curIdialate;
                res.E = max(res.E,E*score);
                if counter>prm.maxNumOfEdges
                    return;
                end
            end
        end
    end
    fprintf('EdgesBeforeNMS = %d\nEdgesAfterNMS = %d\n',length(SC),counter)
    
end

function E = addEdge(ind0,ind1,tSize,E,prm,level)
    global res;
    
    if level == 50
        E = -1;
        disp('deep');
        return;
    end
    
    s0 = res.S0(ind0,ind1);
    res.S0(ind0,ind1) = 0;
    %fprintf('%d,%d,%d\n',ind0,s0,ind1);
    nullVal = 0;
    
    if prm.complexity == 2;
        s1 = res.S1(ind0,ind1);
    else
        s1 = nullVal;
    end
    
    if s1 == ind0 || s1 == ind1
        s1 = nan;
    end
    if s0 == ind0 || s0 == ind1
        s0 = nullVal;
    end
    
    if s1 ~= nullVal && s0 ~= nullVal
        E = addEdge(ind0,s0,tSize,E,prm,level+1);
        if E == -1; return; end;
        E = addEdge(s0,s1,tSize,E,prm,level+1);
        if E == -1; return; end;
        E = addEdge(s1,ind1,tSize,E,prm,level+1);
        if E == -1; return; end;
    elseif s0 ~= nullVal
        E = addEdge(ind0,s0,tSize,E,prm,level+1);
        if E == -1; return; end;
        E = addEdge(s0,ind1,tSize,E,prm,level+1);
        if E == -1; return; end;
    else
        pixels = res.pixels(ind0,ind1,:);
        pixels = pixels(:);
        pixels(pixels == 0) = [];
        E(pixels) = true;
        
        if (numel(pixels) == 0)
            %disp('zero');
            %s = res.S == ind0 | res.S == ind1;
            %imshow(s);
            %fprintf('L = %2.2f, R = %2.2f, C = %2.2f\n',res.L(ind0,ind1),res.R(ind0,ind1),res.C(ind0,ind1));
            E = -1;
        end
    end
end

%% Sub Functions

function [I,S] = subIm(I,S,x0,y0,x1,y1)
    I = I(x0:x1,y0:y1);
    S = S(x0:x1,y0:y1);
end

function createMergeFile(prefix,fileName,prm)
    c = dir(sprintf('Files%d',prm.complexity));
    
    resSingle = [];
    resDouble = [];
    
    for i=1:length(c)
        name = c(i).name;
        
        if strfind(name,prefix) == 1
            load(sprintf('Files%d/%s',prm.complexity,name));
            resSingle = [resSingle;tableSingle];
            resDouble = [resDouble;tableDouble];
        end
    end
    
    s.tableSingle = resSingle;
    s.tableDouble = resDouble;
    
    if ~isempty(c)
        save(fileName, '-struct', 's','-v7.3');
    end
end

function T = threshold(L)
    global res;
    prm = res.prm;
    alpha = 4;
    beta = prm.complexity*2-1;
    w = 2*prm.w;
    T = prm.sigma*sqrt(2*(log(6*res.N)+0*(beta.*L./alpha)*log(2))./(w.*L));
end

function angle = getAngle(ind0,ind1,n)    
    [x0,y0] = ind2sub([n,n],ind0);
    [x1,y1] = ind2sub([n,n],ind1);    
    v1 = x1-x0;
    v2 = y1-y0;
    
    angle = atan(v2./v1);
    angle = 180*angle/pi;
    
    fix = v1<0;
    angle(fix) = angle(fix)+180;
    angle = mod(angle,360);
end

