function [ pitFloorModelS ] = makeAvgPitFloor( trainingMaps, convMap, configS, figHSA)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

figHS = figHSA(1);
%figHSprev = figHSA(2);

[yLen, xLen] = size(trainingMaps{1}.mapZ);
%center = round(/2);
xc = round(xLen/2);
yc = round(yLen/2);

pitGeomS(length(trainingMaps)) = nullPitGeomS();
for i = 1:length(trainingMaps)
    pitGeomS(i) = nullPitGeomS();
end

%floorM{length(trainingMaps)} = [];
floorM = zeros(yLen, xLen, length(trainingMaps));
for i=1:length(trainingMaps)
    mapS = trainingMaps{i};
    
    %[ pitGeomS(i), boundary{i}, borders{i}, boundaryWithoutOutliers{i}, isFit  ] = ...
    [ pitGeomS(i), ~, ~, ~, ~  ] = ...
        fitOneEllipse( mapS, convMap, mapS.label, [], configS, 1, figHS, [] );
    
    ellipseFilter = makeEllipseFilter(mapS, pitGeomS(i), configS.templateMaskInner);
    %ellipseFilter(ellipseFilter == 0) = NaN;
    
    floorMap = mapS.mapZ .* ellipseFilter;
    floorMap(isnan(floorMap)) = 0;
    
    % Need to use round numbers otherwise the image gets "mixed" with the
    % zero-valued cells from the Filter and the nan's
    floorMap = mytranslate(floorMap,round([xc-pitGeomS(i).ellipseXc, yc-pitGeomS(i).ellipseYc]));
    floorMap(floorMap == 0) = nan;
    
    %myImageSc(floorMap);
    
    floorM(:,:,i) = floorMap;
    
end

% I only want to include pixels if several maps include an elevation for
% that pixel. Otherwise outliers could be included in the model.
countMap = ~isnan(floorM);
countMap = sum(countMap, 3);

% Make an average map, but consider only pixels where all maps have data.
% First do a non-weighted average.
averageMap0 = sum(floorM, 3);
averageMap0 = averageMap0 ./ countMap;
countMapTrim = countMap;
countMapTrim(countMap < configS.pitFloorOverlapN) = nan;
averageMap0 = averageMap0 - mean(averageMap0(~isnan(averageMap0)));
%myImageSc(averageMap0);

% Now try to line up each individual floor map with the average
%delta = zeros(length(traingMaps));
for i=1:length(trainingMaps)
    distM = floorM(:,:,i) - averageMap0;
    delta = mean(distM(~isnan(distM)));
    floorM(:,:,i) = floorM(:,:,i) - delta;
end

% Calcuate an average map, but use Thompson's tau test to see if any of the
% elevations are outliers. Get rid of outliers, but only so long as the
% number of points remaining are at most one less than configS.pitFloorOverlapN
% otherwise, just set that pixel to nan.

floorMZero = floorM;
floorMZero(isnan(floorMZero)) = 0;
%averageMap = nan*zeros(size(countMapTrim));
averageMap = sum(floorMZero, 3);
averageMap = averageMap ./ countMapTrim;
%myImageSc(averageMap);

%stdMap = nan*zeros(size(countMapTrim));
stdMap = (floorM - repmat(averageMap,[1,1,size(floorM,3)])).^2;
stdMap(isnan(stdMap)) = 0;
stdMap = sum(stdMap, 3);
stdMap = stdMap ./ (countMapTrim-1);
stdMap = sqrt(stdMap);

%myImageSc(isnan(averageMap) + 2*isnan(stdMap));

%outlierCount = nan*zeros(size(countMapTrim));
%{
[rId, cId] = find( ~isnan(countMapTrim) ) ;
for i=1:length(rId)
    elevs = floorM(rId(i),cId(i),:);
    elevs = elevs(:);
    elevs(isnan(elevs)) = [];
    %indOut = find_outliers_Thompson(elevs, 0.01);
    %indOut = [];
    %elevs(indOut) = [];
    if length(elevs) >= configS.pitFloorOverlapN - 1
        averageMap(rId(i),cId(i)) = mean(elevs);
        stdMap(rId(i),cId(i)) = std(elevs);
    else
        averageMap(rId(i),cId(i)) = nan;
        stdMap(rId(i),cId(i)) = nan;
    end
    outlierCount(rId(i), cId(i)) = length(indOut);
end
%}
%myImageSc(averageMap);
%myImageSc(stdMap);
%myImageSc(outlierCount);

% Outlier detection; get rid of pixels that contribute to a large standard
% deviation.
stdVals = stdMap(~isnan(stdMap));
avgStdMap = mean(stdVals);
sdStdMap = std(stdVals);
averageMap(stdMap > 2*sdStdMap + avgStdMap) = nan;
stdMap(stdMap > 2*sdStdMap+ avgStdMap) = nan;
averageMap = averageMap - mean(averageMap(~isnan(averageMap)));


%disp('If this is anything aside from 0 you need to fix the makeAvgPitFloor code');
%sum(sum(averageMap == 0))



scrsz = get(0,'ScreenSize');
set(0,'Units','pixels');
% Outerposition: lower left corner X,Y, width, height
sumFigH = figure('OuterPosition',...
    [0*scrsz(3) 0.1*scrsz(4) scrsz(3) 0.8*scrsz(4)], 'Visible', 'on');
subplot(1,3,1);
suptitle('Model pit floor');


maxZ = 2*std(averageMap(~isnan(averageMap)));
        
mapH = trainingMaps{1};
mapH.mapZ = averageMap;
pitImageSc(mapH, 2*[-maxZ maxZ]);
hold on;

B = bwboundaries(~isnan(averageMap));

minX = [];
maxX = [];
minY = [];
maxY = [];
for k=1:length(B)
    boundary = B{k};
    if isempty(minX) || min(boundary(:,2)) < minX
        minX = min(boundary(:,2));
    end
    if isempty(maxX) || max(boundary(:,2)) > maxX
        maxX = max(boundary(:,2));
    end
    if isempty(minY) || min(boundary(:,1)) < minY
        minY = min(boundary(:,1));
    end
    if isempty(maxY) || max(boundary(:,1)) > maxY
        maxY = max(boundary(:,1));
    end
    %boundarySize(k) = (max(boundary(:,1)) - min(boundary(:,1))) + ...
    %                  (max(boundary(:,2)) - min(boundary(:,2)));
    plot(boundary(:,2)* mapS.scale, boundary(:,1)* mapS.scale, 'm-');
end


dX = (maxX - minX);
dY = (maxY - minY);
if dX/4 > dY/3
    axisUnit = dX/4;
else
    axisUnit = dY/3;
end
axisUnit = axisUnit*1.2;
xlim([xc - 2*axisUnit xc + 2*axisUnit] * mapS.scale);
ylim([yc - 1.5*axisUnit yc + 1.5*axisUnit] * mapS.scale);
title('Model/average pit floor (nm)');


subplot(1,3,2);

mapH.mapZ = stdMap;
pitImageSc(mapH);%, 0, 'hot');
xlim([xc - 2*axisUnit xc + 2*axisUnit] * mapS.scale);
ylim([yc - 1.5*axisUnit yc + 1.5*axisUnit] * mapS.scale);
title('Model/average pit floor standard deviation (nm)');

subplot(1,3,3); 
mapH.mapZ = stdMap;
hist(stdVals, 200);
title('Standard deviation histogram (nm)');


button = questdlg('Do you like the look of the model pit floor?', 'Yes', 'No');
if ~strcmp(button, 'Yes')
    error('You said No!');
end

close(sumFigH);


pitFloorModelS = struct(...
    'model', averageMap,...
    'modelSD', stdMap,...
    'minX', minX, ...
    'maxX', maxX, ...
    'minY', minY, ...
    'maxY', maxY, ...
    'dX', dX, ...
    'dY', dY);


%myImageSc(xor(isnan(pitFloorModelS.model), isnan(pitFloorModelS.modelSD)));
%error('It is possible that some of the model pixels can be real numbers while the corresponding SD map has NANs');

pitFloorModelS.boundaries = B;

end
