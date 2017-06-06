function [] = plotPatchesHeights(W2d,points,ColourMatrix,npoints)
% Checks which points overlap and plots them over W2d.

Nm = length(points);

if ~exist('ColourMatrix','var')
    ColourMatrix = rand(Nm,3);
else
    if isempty(ColourMatrix)
        ColourMatrix = rand(Nm,3);
    end
end


if ~exist('npoints','var')
    npoints = [];
end

dh = 5*(max(W2d(:)) - min(W2d(:)))/Nm;
h = zeros(1,size(W2d,2));


plot3(W2d(1,points{1}),W2d(2,points{1}),h(1,points{1}),'ko','MarkerSize',10,'MarkerFaceColor',ColourMatrix(1,:),'Color',ColourMatrix(1,:),'LineWidth',1)
hold on
for m = 2:Nm
    h = h + dh;
    p = points{m};
    if ~isempty(p)
        plot3(W2d(1,p),W2d(2,p),h(1,p),'ko','MarkerSize',10,'MarkerFaceColor',ColourMatrix(m,:),'Color',ColourMatrix(m,:),'LineWidth',1)
    end
    %     m
    %     pause();
end

if ~isempty(npoints)
    h =  h + dh;
    plot3(W2d(1,npoints),W2d(2,npoints),h(1,npoints),'k+','MarkerSize',10,'LineWidth',1)
end
axis equal
hold off