% -----------------------------------------
% Graph-Cuts for F-Formation (GCFF)
% 2015 - University of Verona
% Written by Francesco Setti
% -----------------------------------------
%
% This function deletes all the groups formed by a single individual.
%

function gr_out = ff_deletesingletons( gr_in )


% Delete all the groups in gr_in with a single component.
tokeep = [] ;
for idxGroup = 1:length(gr_in)
    if length(gr_in{idxGroup})>1
        tokeep = [tokeep, idxGroup] ;
    end
end
gr_out = gr_in(tokeep) ;
