function [ elevDeepest, elevDeepestSD ] = pitFloorElevation( mapS, configS, pitGeomS, pfmS, figHS )
%PITELEVATION Summary of this function goes here
%   Detailed explanation goes here

%  ** Get the average elevation in the centre of the pit
%  - User input (in config file) for diameter
%  - Get average elevation in circle
%pitFloorMeasRadius = configS.pitFloorMeasuredDiameter / mapS.scale / 2;
 
%[colsInImage, rowsInImage] = meshgrid(1:length(Map(1,:)), 1:length(Map(:,1)));
%circleMask = (mapS.y - yc).^2 + (mapS.x - xc).^2 < (pitFloorMeasRadius)^2;



center = round(size(mapS.mapZ)/2);
xc = center(2);
yc = center(1);
tx = round(xc-pitGeomS.ellipseXc);
ty = round(yc-pitGeomS.ellipseYc);

if isnan(tx) || isnan(ty)
    disp('tx or ty are NaN in pitFloorElevation');
    elevDeepest = nan;
    elevDeepestSD = nan;
    return;
end
translatedMap = mytranslate(mapS.mapZ,[tx ty]);
%myImageSc(translatedMap);

maskedMap = translatedMap .* ~isnan(pfmS.model);
maskedMap(maskedMap == 0) = nan;

%{
if sum(~isnan(maskedMap(:))) == 0
    disp('The modelPitFloor and the actual pit floor do not overlap at all!');
    elevDeepest = nan;
    elevDeepestSD = nan;
    return;
end
%}
%myImageSc(maskedMap);

bigBlobMap = bwareaopen(~isnan(maskedMap), ...
    round((configS.minimumFloorBlobSize/100)*pitGeomS.ellipseArea/(mapS.scale)^2));
%bigBlobMap(bigBlobMap == 0) = nan;
%myImageSc(bigBlobMap);
cc = bwconncomp(bigBlobMap);
if cc.NumObjects >= 1
    kLen(cc.NumObjects) = 0;
    for i=1:cc.NumObjects
        kLen(i) = length(cc.PixelIdxList{i});
    end
    maxLen = max(kLen);
    % Blobs must be at least x% of the size of the biggest blob to be
    % considered.

    elev = ones(1,cc.NumObjects) * nan;
    elevSD = ones(1,cc.NumObjects) * nan;
    for i=1:cc.NumObjects
        if kLen(i) >= (configS.minimumBlobRelativeArea/100)*maxLen
            blob = zeros(size(bigBlobMap));
            blob(cc.PixelIdxList{i}) = 1;
            blob(blob == 0) = nan;
            %myImageSc(blob);
            blob = blob .* maskedMap;

            distM = blob - pfmS.model;
            %myImageSc(distM);
            distM(isnan(distM)) = 0;
            %pfmS.modelSD(pfmS.modelSD == 0) = nan;
            weights = 1./(pfmS.modelSD).^2;
            weights(isnan(weights)) = 0;
            weights(isnan(blob)) = 0;
            weightedDist = distM .* weights;
            %myImageSc(weights);
            %myImageSc(weightedDist);
            weightedDistAvg = sum(weightedDist(weightedDist ~= 0));
            weightedDistAvg = weightedDistAvg / sum(weights(weightedDist ~= 0));

            elev(i) = weightedDistAvg;
            %wDist = weightedDist / sum(weights(weightedDist ~= 0));
            %myImageSc( wDist );
            %myImageSc(blob - (modelPitFloor + elev(i)));
            % For the standard deviation, see:
            % https://en.wikipedia.org/wiki/Weighted_arithmetic_mean#Dealing_with_variance
            %delta = distM - weightedDistAvg;

            delta = distM - weightedDistAvg;
            delta(distM == 0) = 0;
            if sum(sum(xor((delta == 0), (weights == 0)))) && sum(sum(distM ~= 0)) > 1
                error('There is a tricky problem in pitFloorElevation');
            end
            elevSD(i) = sqrt( sum(sum(weights.*(delta).^2)) / sum(weights(weightedDist ~= 0)) );
        end
    end

    % TODO: Now choose the keepMBlobs deepest ones.
    [elevDeepest,ind] = min(elev);
    elevDeepestSD = elevSD(ind);
else
    elevDeepest = nan;
    elevDeepestSD = nan;
    ind = nan;
end
    



%{
elev = mapS.mapZ .* circleMask;

elev = elev(circleMask ~= 0);
elev = elev(~isnan(elev));
medianElev = median(elev);
sdElev = std(elev);
%}

if isempty(elevDeepest)
    elevDeepest = nan;
    elevDeepestSD = nan;
end

if ~isempty(figHS.summary)
    set(0, 'CurrentFigure', figHS.summary);
    subplot(3,2,5,'replace');
    
    distM = [];
    maxZ = [];
    B1 = [];
    B2 = [];
    if ~isnan(ind)
        blob = zeros(size(bigBlobMap));
        blob(cc.PixelIdxList{ind}) = 1;
        blob(blob == 0) = nan;
        blob = blob .* maskedMap;
        distM = blob - elevDeepest - pfmS.model;
        
        
        maxZ = 2*std(distM(~isnan(distM)));
        
        B1 = bwboundaries(bigBlobMap);
        B2 = bwboundaries(~isnan(blob));
    else
        distM = zeros(size(pfmS.model));
    end
    
    if ~isnan(ind)
        pitImageSc(mapS, elevDeepest + 2*[-maxZ maxZ]);
    else
        pitImageSc(mapS);
    end
    if ~isnan(xc) && ~isnan(yc)
        if pfmS.dX/4 > pfmS.dY/3
            axisUnit = pfmS.dX/4;
        else
            axisUnit = pfmS.dY/3;
        end
        axisUnit = axisUnit*1.2;
        xlim([xc-tx - 2*axisUnit xc-tx + 2*axisUnit] * mapS.scale);
        ylim([yc-ty - 1.5*axisUnit yc-ty + 1.5*axisUnit] * mapS.scale);
    end
    
    hold on;
    

    for k=1:length(pfmS.boundaries)
        boundary = pfmS.boundaries{k};
        plot((boundary(:,2)-tx)*mapS.scale, (boundary(:,1)-ty)* mapS.scale,'m-', 'LineWidth',1);
    end
    for k=1:length(B1)
        boundary = B1{k};
        plot((boundary(:,2)-tx)*mapS.scale, (boundary(:,1)-ty)*mapS.scale, 'g-');
    end
    for k=1:length(B2)
        boundary = B2{k};
        plot((boundary(:,2)-tx)*mapS.scale, (boundary(:,1)-ty)*mapS.scale, 'y-', 'LineWidth',2);
    end
   
    title('Pit floor (nm) with outlines of the blobs and the model floor');

    subplot(3,2,6,'replace');
    %dElev = 0.01E4;
    
    mapH = mapS;
    %mapH.mapZ = mytranslate(distM,-1*[tx ty]);
    mapH.mapZ = mytranslate(translatedMap - elevDeepest - pfmS.model,-1*[tx ty]);
    

    if ~isnan(ind)
        %hist(hDist, min(hDist)-dElev/2:dElev:max(hDist)+dElev/2);%, 'FaceColor', [1 0 0])
        pitImageSc(mapH, [-maxZ maxZ]);
        xlim([xc-tx - 2*axisUnit xc-tx + 2*axisUnit] * mapS.scale);
        ylim([yc-ty - 1.5*axisUnit yc-ty + 1.5*axisUnit] * mapS.scale);
        
        
        %histogram(elev); , 'FaceColor', [1 0 0]
        hold on;
        %vline(h, medianElev, 'r-');
        title('Distance from the pit floor to the fitted model floor (nm)');
        for k=1:length(pfmS.boundaries)
            boundary = pfmS.boundaries{k};
            plot((boundary(:,2)-tx)*mapS.scale, (boundary(:,1)-ty)* mapS.scale,'m-', 'LineWidth',1);
        end
        for k=1:length(B1)
            boundary = B1{k};
            plot((boundary(:,2)-tx)*mapS.scale, (boundary(:,1)-ty)*mapS.scale, 'g-');
        end
        for k=1:length(B2)
            boundary = B2{k};
            plot((boundary(:,2)-tx)*mapS.scale, (boundary(:,1)-ty)*mapS.scale, 'y-', 'LineWidth',2);
        end
    end
    %}
end

end

