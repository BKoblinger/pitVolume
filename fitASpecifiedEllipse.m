function [ xc, yc, boundaryTrim ] = fitASpecifiedEllipse( mapS, xc0, yc0, a, b, phi, convMap, configS )
%FITASPECIFIEDELLIPSE Summary of this function goes here
%   Detailed explanation goes here
approxPitRadiusPixels = configS.pitDiameter / mapS.scale / 2;
approxPitFloorRadiusPixels = configS.pitFloorDiameter / mapS.scale / 2;

ellipseFilter = ...
        (((mapS.x-(max(max(mapS.x)/2)))*cos(phi)+(mapS.y-(max(max(mapS.y)/2)))*sin(phi)).^2/(a/mapS.scale)^2) + ...
        (((mapS.x-(max(max(mapS.x)/2)))*sin(phi)-(mapS.y-(max(max(mapS.y)/2)))*cos(phi)).^2/(b/mapS.scale)^2) <= 1;
    
ellipseFilterEdge = edge(ellipseFilter);
ellipseFilterEdge = ellipseFilterEdge / sum(sum(ellipseFilterEdge));
gaussFilter = fspecial('gaussian', 6, 1.5);
if configS.showImgAnalFigs == 1
    myImageSc(ellipseFilter);
    myImageSc(ellipseFilterEdge);
    myImageSc(gaussFilter);
end

ellipseFilterEdge = imfilter(ellipseFilterEdge, gaussFilter);
ellipseFilterEdge( ~any(ellipseFilterEdge,2), : ) = [];  %rows
ellipseFilterEdge( :, ~any(ellipseFilterEdge,1) ) = [];  %columns

if configS.showImgAnalFigs == 1
    myImageSc(ellipseFilterEdge);
end


% Set to zero anything from the center of the previously found ellipse
%[filterSize] = size(ellipseFilter);

%nanMap(1:(yc0-filterSize(1)),:) = 0;
%nanMap(yc0+filterSize(1):end,:) = 0;
%myImageSc(nanMap);
%nanMap(:,1:(xc0-filterSize(2))) = 0;
%nanMap(:,xc0+filterSize(2):end) = 0;
%myImageSc(nanMap);

%borderM = zeros(size(mapS.mapZ));
%borderM(sub2ind(size(borderM),boundary(:,1),boundary(:,2))) = 1;
borderM = isnan(mapS.mapZ);
if ~isempty(convMap)
    [ edgeMap ] = deNoise( borderM, mapS.scale, convMap, 0, 0,  configS );
%edgeMap = deNoise(borderM, approxPitRadiusPixels, ...
%    approxPitFloorRadiusPixels, configS.showImgAnalFigs);
else
    edgeMap = edge(borderM);
end
corrM = imfilter(1.0*edgeMap, ellipseFilterEdge);

if configS.showImgAnalFigs == 1
    myImageSc(borderM);
    myImageSc(edgeMap);
    myImageSc(corrM);
end

[yc, xc] = find(corrM == max(max(corrM))); 

ellipseFilter = ...
        (((mapS.x-(xc))*cos(phi)+(mapS.y-(yc))*sin(phi)).^2/(a/mapS.scale)^2) + ...
        (((mapS.x-(xc))*sin(phi)-(mapS.y-(yc))*cos(phi)).^2/(b/mapS.scale)^2) <= 1;
ellipseFilterEdge = 1.0*edge(ellipseFilter);
widenFilter = fspecial('disk', 5);
if configS.showImgAnalFigs == 1
    myImageSc(widenFilter);
end
fatEllipse = imfilter(ellipseFilterEdge, widenFilter);
fatEllipse(fatEllipse ~= 0) = 1;
if configS.showImgAnalFigs == 1
    myImageSc(fatEllipse);
end

boundaryTrimM = edgeMap .* fatEllipse;
if configS.showImgAnalFigs == 1
    myImageSc(boundaryTrimM);
end

[x, y] = find(boundaryTrimM);
boundaryTrim = [x y];


end

