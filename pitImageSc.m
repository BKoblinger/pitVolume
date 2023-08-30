function [  ] = pitImageSc( mapS, lim, cmS )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

cm = flipud(brewermap(64,'RdBu'));
if nargin == 3 && strcmp(cmS, 'hot')
    cm = hot;
end

if nargin == 1 || (nargin == 3 && lim == 0)
    imagesc(mapS.x(1,:) * mapS.scale, mapS.y(:,1) * mapS.scale, mapS.mapZ);
elseif nargin == 2
    if any(isnan(lim)) || any(isinf(lim)) || lim(1) == lim(2)
        imagesc(mapS.x(1,:) * mapS.scale, mapS.y(:,1) * mapS.scale, mapS.mapZ);
    else
        imagesc(mapS.x(1,:) * mapS.scale, mapS.y(:,1) * mapS.scale, mapS.mapZ, lim);
    end
elseif nargin >= 4
    error('Invalid number of arguments passed to pitImageSc');
end
set(gca, 'YDir', 'normal');
axis equal;
axis tight;
colormap(cm);
colorbar;


end

