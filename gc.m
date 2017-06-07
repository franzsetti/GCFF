function seg = gc( f, stride, MDL_in)
% Runs graph-cuts


locs = find_locs( f, stride ) ;
unary = init( locs, f, MDL_in ) ;
neigh  = zeros(size(f,1),1) ;
weight = zeros(size(f,1),1) ;
seg = (0:size(f,1)-1)' ;

segold = zeros(length(seg),1) ;

MAX_ITER = 10 ;
numiter = 1 ;

% close all,
while ~isequal(seg,segold) && numiter <= MAX_ITER
%     %%%%%% TEMP
%     figure(numiter),
%     ff_plot_person_tv(f,'b',50),
%     hold on,
%     plot(locs(:,1),locs(:,2),'b*')
%     for ii = 1:size(locs,1)
%         text(locs(ii,1),locs(ii,2),['C_{',num2str(ii),'}'], ...
%             'Color','b','FontWeight','bold','HorizontalAlignment','center') ;
%     end
%     labs = seg ;
%     ulabs = unique(labs) ;
%     clear plocs;
%     for ii = 1:length(ulabs)
%         plocs(ii,:) = mean(locs(labs==ulabs(ii),:), 1) ;
%     end
%     plot(plocs(:,1),plocs(:,2),'r*')
%     for ii = 1:size(plocs,1)
%         text(plocs(ii,1),plocs(ii,2),['O_{',num2str(ii),'}'], ...
%             'Color','r','FontWeight','bold','HorizontalAlignment','center') ;
%     end
%     axis equal
% %     ylim([-400,250])
% %     xlim([-300,300])
% %     if numiter>2
%         keyboard
% %     end
%     %%%%%% END TEMP

    segold = seg ;
    mdl = ones(size(unary,2),1) * MDL_in ;
    % run graphcuts
    [aa,seg] = expand( unary, neigh', weight', seg', mdl, Inf ) ;
    % discard unused labels
    %     seg = uint8(seg) ;
    %     seg = unique( seg ) ;
    % refit distances
    [unary,seg] = calc_distance( locs, f, seg, MDL_in) ;
    numiter = numiter + 1 ;
end

end  % function gc


% -----------------------
%  Additional functions
% -----------------------

function locs = find_locs( f, stride )
% Estimate focal centers for each person given features

locs = zeros( size(f,1), 2) ;
locs(:,1) = f(:,2) + cos(f(:,4))*stride ;
locs(:,2) = f(:,3) + sin(f(:,4))*stride ;

end  % function find_locs


function distmat = init( locs, f, mdl )

distmat = calc_distance( locs, f, 1:size(locs,1), mdl ) ;

end  % function init


function [distmat,labels] = calc_distance( loc, f, labels, mdl)
% Given focal locations, raw locations(f) and initial labelling l find cost
% of assigning people to new locations given by the mean of their
% labelling.

u = unique( labels ) ;
distmat = zeros( size(loc,1), length(u) ) ;

for i = 1:length(u)

    means = mean( loc(labels==u(i),:), 1) ;
    labels(labels==u(i)) = i-1 ;
    % dist(:,i)=((loc-means)**2).sum(2) %%you need to use repmat here, or some kind of built in L_2 norm on a matrix of vectors
    distmat(:,i) = pdist2(loc,means).^2 ;
    %computed sum-squares distance, now
    mask = find( distmat(:,i)<mdl ) ;
    disp = f(:,2:3) - repmat(means,[size(distmat,1),1]) ;
    for j = mask'
        for k = mask'
            distk = norm(disp(k,:)) ;
            distj = norm(disp(j,:)) ;
            if distk > distj
                inner = disp(k,:) * disp(j,:)' ;
                norma = distk * distj ;
                if inner/norma > 0.75
                    distmat(k,i) = distmat(k,i) + 100^(inner/norma*distk/distj) ;
                end
            end
        end
    end
end

end  % function calc_distance
