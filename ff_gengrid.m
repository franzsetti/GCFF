% -----------------------------------------
% Graph-Cuts for F-Formation (GCFF)
% 2015 - University of Verona
% Written by Francesco Setti
% -----------------------------------------
%
% This function generates a grid for Hough voting with a quantization
% 'quant'.
%

function [votegrid, votegrid_pos] = ff_gengrid( features, param, quant)


if ~exist('quant','var')
    quant = 1 ;
end

% Initialize boundaries
idxFrame = 1 ;
while isempty(features{idxFrame})
    idxFrame = idxFrame + 1 ;
end
x_min = min(features{idxFrame}(:,2)) ;
x_max = max(features{idxFrame}(:,2)) ;
y_min = min(features{idxFrame}(:,3)) ;
y_max = max(features{idxFrame}(:,3)) ;

% Update boundaries
for ii = 1:length(features)
    if ~isempty(features{ii})
        x_min = min(x_min,min(features{ii}(:,2))) ;
        x_max = max(x_max,max(features{ii}(:,2))) ;
        y_min = min(y_min,min(features{ii}(:,3))) ;
        y_max = max(y_max,max(features{ii}(:,3))) ;
    end
end

% Set boundaries of grid
xmin = x_min - 2 * param.radius ;
xmax = x_max + 2 * param.radius ;
ymin = y_min - 2 * param.radius ;
ymax = y_max + 2 * param.radius ;

% Generate grid
[votegrid(:,:,1),votegrid(:,:,2)] = meshgrid(xmin:quant:xmax,ymin:quant:ymax) ;

% Generate votegrid_pos
votegrid_pos = reshape(votegrid,[],2) ;
