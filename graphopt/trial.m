load Kmatrix.mat
 [a b]=allgc(K,s',zeros(4800,1),zeros(4800,1),ones(4800,1),1, ...
             zeros(165,1),Inf);
 
  [a b]=allgc(K,s',zeros(4800,1),zeros(4800,1),ones(4800,1),1,45* ...
              ones(165,1),Inf);
  
    [a b]=allgc(K,s',s',ones(size(s')),ones(4800,1),1,45*ones(165, ...
                                                      1),Inf);
    
    [a b]=allgc(K,s',s',ones(size(s')),ones(4800,1),1,450*ones(165, ...
                                                      1),Inf);
    
    
