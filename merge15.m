function [tableSingleFinal,tableDoubleFinal] = merge15(S11,S21,S12,S22)
    n = size(S11,1);
    disp(n);
    tableSingle = nan(40*n^3,3);
    indexSingle = 1;
    tableDouble = nan(40*n^3,3);
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
    
    x = n; y = 1:n; split5  = [S11(x,y) S21(x,y(2:end))];
    x = 1; y = 1:n; up  = [S11(x,y) S21(x,y(2:end))];
    x = n; y = 1:n; down  = [S12(x,y) S22(x,y(2:end))];
    
    x = 1:n; y = n; split6  = [S11(x,y)' S12(x(2:end),y)'];
    x = 1:n; y = 1; left  = [S11(x,y)' S12(x(2:end),y)'];
    x = 1:n; y = n; right  = [S21(x,y)' S22(x(2:end),y)'];
    
    disp('Upper Rectangle')

    disp('e1');   
    [tableSingle,indexSingle] = addIndices(tableSingle,indexSingle,e1,split1,e3);
    [tableSingle,indexSingle] = addIndices(tableSingle,indexSingle,e1,split1,e4);
    [tableSingle,indexSingle] = addIndices(tableSingle,indexSingle,e1,split1,split2);
    
    disp('e2');
    [tableSingle,indexSingle] = addIndices(tableSingle,indexSingle,e2,split1,e4);
    [tableSingle,indexSingle] = addIndices(tableSingle,indexSingle,e2,split1,split2);
    
    disp('s4');
    [tableSingle,indexSingle] = addIndices(tableSingle,indexSingle,split4,split1,e3);
    [tableSingle,indexSingle] = addIndices(tableSingle,indexSingle,split4,split1,e4);
        
    disp('Lower Rectangle')
    
    disp('e8');
    [tableSingle,indexSingle] = addIndices(tableSingle,indexSingle,e8,split3,e6);
    [tableSingle,indexSingle] = addIndices(tableSingle,indexSingle,e8,split3,e5);
    [tableSingle,indexSingle] = addIndices(tableSingle,indexSingle,e8,split3,split2);
    
    disp('e7');
    [tableSingle,indexSingle] = addIndices(tableSingle,indexSingle,e7,split3,e5);
    [tableSingle,indexSingle] = addIndices(tableSingle,indexSingle,e7,split3,split2);
    
    disp('s4');
    [tableSingle,indexSingle] = addIndices(tableSingle,indexSingle,split4,split3,e6);
    [tableSingle,indexSingle] = addIndices(tableSingle,indexSingle,split4,split3,e5);
    
    %disp('Left Rectangle')

    %disp('e2');   
    %[tableSingle,indexSingle] = addIndices(tableSingle,indexSingle,e2,split4,e7);
    %[tableSingle,indexSingle] = addIndices(tableSingle,indexSingle,e2,split4,e8);
    %[tableSingle,indexSingle] = addIndices(tableSingle,indexSingle,e2,split4,split3);
    
    %disp('e1');
    %[tableSingle,indexSingle] = addIndices(tableSingle,indexSingle,e1,split4,e7);
    %[tableSingle,indexSingle] = addIndices(tableSingle,indexSingle,e1,split4,split3);
    
    %disp('s4');
    %[tableSingle,indexSingle] = addIndices(tableSingle,indexSingle,split1,split4,e7);
    %[tableSingle,indexSingle] = addIndices(tableSingle,indexSingle,split1,split4,e8);
        
    %disp('Right Rectangle')
    
    %disp('e3');
    %[tableSingle,indexSingle] = addIndices(tableSingle,indexSingle,e3,split2,e5);
    %[tableSingle,indexSingle] = addIndices(tableSingle,indexSingle,e3,split2,e6);
    %[tableSingle,indexSingle] = addIndices(tableSingle,indexSingle,e3,split2,split3);
    
    %disp('e4');
    %[tableSingle,indexSingle] = addIndices(tableSingle,indexSingle,e4,split2,e6);
    %[tableSingle,indexSingle] = addIndices(tableSingle,indexSingle,e5,split2,split3);
    
    %disp('s1');
    %[tableSingle,indexSingle] = addIndices(tableSingle,indexSingle,split1,split2,e5);
    %[tableSingle,indexSingle] = addIndices(tableSingle,indexSingle,split1,split2,e6);
    
    disp('Whole Square')
 
    %disp('left');
    %[tableDouble,indexDouble] = addIndices(tableDouble,indexDouble,left,split6,e3);
    %[tableDouble,indexDouble] = addIndices(tableDouble,indexDouble,left,split6,right);
    %[tableDouble,indexDouble] = addIndices(tableDouble,indexDouble,left,split6,e8);

    disp('up');
    [tableDouble,indexDouble] = addIndices(tableDouble,indexDouble,up,split5,e5);
    [tableDouble,indexDouble] = addIndices(tableDouble,indexDouble,up,split5,down);
    [tableDouble,indexDouble] = addIndices(tableDouble,indexDouble,up,split5,e8);
    
    disp('e1');
    [tableDouble,indexDouble] = addIndices(tableDouble,indexDouble,e1,split5,e5);
    [tableDouble,indexDouble] = addIndices(tableDouble,indexDouble,e1,split5,down);
    
    disp('e4');
    [tableDouble,indexDouble] = addIndices(tableDouble,indexDouble,e4,split5,down);
    [tableDouble,indexDouble] = addIndices(tableDouble,indexDouble,up,split5,e8);
    
    tableSingleFinal = tableSingle(1:indexSingle-1,:);
    tableDoubleFinal = tableDouble(1:indexDouble-1,:);
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
