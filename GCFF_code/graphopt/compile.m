% mex  -v CXXFLAGS="-O3 -ffast-math  \$CXXFLAGS" allgc.cpp  
% mex  -v CXXFLAGS="-O3 -ffast-math  \$CXXFLAGS" expand.cpp 
% mex  -v CXXFLAGS="-O3 -ffast-math  \$CXXFLAGS" multi.cpp 
%mex   -v CXXFLAGS="-O3 -ffast-math  \$CXXFLAGS" convex.cpp 
mex allgc.cpp 
mex expand.cpp 
mex multi.cpp 


%load("/homes/jrkf/Public/pre_gc.mat", "K","s","old_labels","Lambda","Np","MDL")