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

function [votegrid, votegrid_pos] = ff_gengrid( features, param, quant)
% This function generates a grid for Hough voting with a quantization
% 'quant'.


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

