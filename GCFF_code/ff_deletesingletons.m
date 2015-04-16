% Copyright 2014 Francesco Setti and Marco Cristani
% Department of Computer Science
% University of Verona
% http://www.di.univr.it/
% 
% Permission is hereby granted, free of charge, to any person obtaining a
% copy of this software and associated documentation files (the
% "Software"), to deal in the Software without restriction, subject to the
% following conditions: 
%  * The above copyright notice and this permission notice shall be
%    included in all copies or substantial portions of the Software. 
%  * The Software is provided "as is", without warranty of any kind.
%
% September 2014
%

function gr_out = ff_deletesingletons( gr_in )
% This function deletes all the groups formed by a single individual.

% Delete all the groups in gr_in with a single component.
tokeep = [] ;
for idxGroup = 1:length(gr_in)
    if length(gr_in{idxGroup})>1
        tokeep = [tokeep, idxGroup] ;
    end
end
gr_out = gr_in(tokeep) ;