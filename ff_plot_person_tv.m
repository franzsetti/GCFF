% -----------------------------------------
% Graph-Cuts for F-Formation (GCFF)
% 2015 - University of Verona
% Written by Francesco Setti
% -----------------------------------------
%
% This function plots a person symbol from top view.
% 

function ff_plot_person_tv( xx_feat, yy_feat, alpha_feat, col, scale)


if ~exist('col','var')
    col = 'b' ;
end

if ~exist('scale','var')
    scale = 1 ;
end

if size(xx_feat,2) == 4
    xx = xx_feat(:,2) ;
    yy = xx_feat(:,3) ;
    alpha = xx_feat(:,4) ;
    ID = xx_feat(:,1) ;
    if nargin > 1
        col = yy_feat ;
    end
    if nargin > 2
        scale = alpha_feat ;
    end
elseif min(size(xx_feat)) == 1
    xx = xx_feat ;
    yy = yy_feat ;
    alpha = alpha_feat ;
end


% Define width and height of person from top view
bw = 0.5 * scale ;   % [m]
bh = 0.2 * scale ;   % [m]
hd = 0.15 * scale ;  % [m]

hold on,
for ii = 1:length(xx)
    ellipse( bw, bh, alpha(ii)+pi/2, xx(ii), yy(ii), col) ;
    ellipse( hd, hd, alpha(ii)+pi/2, xx(ii), yy(ii), col) ;
    %arrow([xx(ii),yy(ii)],[xx(ii)+bw*cos(alpha(ii)),yy(ii)+bw*sin(alpha(ii))]) ;
    quiver(xx(ii),yy(ii),bw*cos(alpha(ii)),bw*sin(alpha(ii)),col)
    if exist('ID')
        % text(xx(ii)+0.1*scale,yy(ii)+0.1*scale,num2str(ID(ii)),'Color',col) ;
        text(xx(ii)-(bh*1.5)*cos(alpha(ii)),yy(ii)-(bh*1.5)*sin(alpha(ii)),['P_{',num2str(ID(ii)),'}'],...
            'Color',col,'FontWeight','bold','HorizontalAlignment','center') ;
    end
end
hold off
