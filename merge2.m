function [tableSingleFinal,tableDoubleFinal1,tableDoubleFinal2] = merge2(S11,S21,S12,S22)
    n = size(S11,1);
    
    tableSingle = zeros(20*n^3,3);
    indexSingle = 1;
    tableDouble = zeros(15*n^4,4);
    indexDouble = 1;
    
    % S = [S11 S21; S12 S22];
    
    % indices initialization
    x = 1:n; y = 1; e1 = S11(x,y);
    x = 1; y = 1:n; e2 = S11(x,y);
    x = 1; y = 1:n; e3 = S21(x,y);
    x = 1:n; y = n; e4 = S21(x,y);
    x = 1:n; y = n; split1 = S11(x,y);
    x = n; y = 1:n; split2 = S21(x,y);
    x = 1:n; y = n; split3 = S12(x,y);
    x = n; y = 1:n; split4 = S11(x,y);
    x = 1:n; y = n; e5 = S22(x,y);
    x = n; y = 1:n; e6 = S22(x,y);
    x = n; y = 1:n; e7 = S12(x,y);
    x = 1:n; y = 1; e8 = S12(x,y);
    
    disp('e1');   
    [tableSingle,indexSingle] = addIndices(tableSingle,indexSingle,e1,split1,e3);
    [tableSingle,indexSingle] = addIndices(tableSingle,indexSingle,e1,split1,e4);
    [tableSingle,indexSingle] = addIndices(tableSingle,indexSingle,e1,split4,e7);
        
    [tableDouble,indexDouble] = addIndices2(tableDouble,indexDouble,e1,split1,split2,e5);
    [tableDouble,indexDouble] = addIndices2(tableDouble,indexDouble,e1,split1,split2,e6);
    [tableDouble,indexDouble] = addIndices2(tableDouble,indexDouble,e1,split4,split3,e5);
    [tableDouble,indexDouble] = addIndices2(tableDouble,indexDouble,e1,split4,split3,e6);
   
    disp('e2');
    [tableSingle,indexSingle] = addIndices(tableSingle,indexSingle,e2,split1,e4);
    [tableSingle,indexSingle] = addIndices(tableSingle,indexSingle,e2,split4,e7);
    [tableSingle,indexSingle] = addIndices(tableSingle,indexSingle,e2,split4,e8);

    [tableDouble,indexDouble] = addIndices2(tableDouble,indexDouble,e2,split1,split2,e5);
    [tableDouble,indexDouble] = addIndices2(tableDouble,indexDouble,e2,split1,split2,e6);
    [tableDouble,indexDouble] = addIndices2(tableDouble,indexDouble,e2,split4,split3,e5);
    [tableDouble,indexDouble] = addIndices2(tableDouble,indexDouble,e2,split4,split3,e6);
    
    disp('e3');
    [tableSingle,indexSingle] = addIndices(tableSingle,indexSingle,e3,split2,e5);
    [tableSingle,indexSingle] = addIndices(tableSingle,indexSingle,e3,split2,e6);

    [tableDouble,indexDouble] = addIndices2(tableDouble,indexDouble,e3,split2,split3,e7);
    [tableDouble,indexDouble] = addIndices2(tableDouble,indexDouble,e3,split2,split3,e8);
    [tableDouble,indexDouble] = addIndices2(tableDouble,indexDouble,e3,split1,split4,e7);
    [tableDouble,indexDouble] = addIndices2(tableDouble,indexDouble,e3,split1,split4,e8);
    
    tableDouble2 = tableDouble(1:indexDouble-1,:);
    tableDouble = zeros(15*n^4,4);
    indexDouble = 1;
    
    disp('e4');
    [tableSingle,indexSingle] = addIndices(tableSingle,indexSingle,e4,split2,e6);

    [tableDouble,indexDouble] = addIndices2(tableDouble,indexDouble,e4,split2,split3,e7);
    [tableDouble,indexDouble] = addIndices2(tableDouble,indexDouble,e4,split2,split3,e8);
    [tableDouble,indexDouble] = addIndices2(tableDouble,indexDouble,e4,split1,split4,e7);
    [tableDouble,indexDouble] = addIndices2(tableDouble,indexDouble,e4,split1,split4,e8);
    
    disp('e5');
    [tableSingle,indexSingle] = addIndices(tableSingle,indexSingle,e5,split3,e7);
    [tableSingle,indexSingle] = addIndices(tableSingle,indexSingle,e5,split3,e8);
    
    disp('e6');
    [tableSingle,indexSingle] = addIndices(tableSingle,indexSingle,e6,split3,e8);
    
    tableSingleFinal = tableSingle(1:indexSingle-1,:);
    tableDoubleFinal1 = tableDouble2;
    tableDoubleFinal2 = tableDouble(1:indexDouble-1,:);
end

function [table,index] = addIndices(table,index,ind0,s0,ind1)
    [I0,S0,I1] = ndgrid(ind0,s0,ind1);
    I0 = I0(:);I1 = I1(:);S0 = S0(:);
    Z = [I0 I1 S0];
    bad = S0 == I0 | S0 == I1 | I0 == I1;
    Z(bad,:) = [];
    [m,n] = size(Z);
    table(index:index+m-1,:) = Z;
    index = index+m;
end

function [table,index] = addIndices2(table,index,ind0,s0,s1,ind1)
    [I0,S0,S1,I1] = ndgrid(ind0,s0,s1,ind1);
    I0 = I0(:);I1 = I1(:);S0 = S0(:); S1 = S1(:);
    Z = [I0 I1 S0 S1];
    bad = I0 == S0 | I0 == S1 | I0 == I1 | S0 == S1 | S0 == I1 | S1 == I1;
    Z(bad,:) = [];
    [m,n] = size(Z);
    table(index:index+m-1,:) = Z;
    index = index+m;
end