function [ edgeMap, mask ] = deNoise( nanMap, scale, templateMap, approxOuterRadiusPixels, approxInnerRadiusPixels,  configS )
%DENOISE Summary of this function goes here
%   Detailed explanation goes here

mask = [];

if configS.showImgAnalFigs
    myImageSc(nanMap);
end

if isempty(templateMap)
    minBlobSize = pi*(approxOuterRadiusPixels)^2 - pi*(approxInnerRadiusPixels)^2;
    minBlobSize = minBlobSize/1.5;
    nanMap = bwareaopen(nanMap, round(minBlobSize));
    nanMapClear = imclearborder(nanMap);
    if configS.showImgAnalFigs
        myImageSc(nanMapClear);
    end

    if sum(nanMapClear(:)) == 0
        nanMapClear = nanMap;
    end

    [labelMap, numObjects] = bwlabel(nanMapClear, 8);

    if configS.showImgAnalFigs
        myImageSc(labelMap);
    end

    stats = regionprops(labelMap,nanMap,'Centroid');

    distToCenter = zeros(1,numObjects);
    for i=1:numObjects
        distToCenter(i) = sum((stats(i).Centroid - 0.5*fliplr(size(nanMap))).^2);
    end

    [~,ind] = min(distToCenter);

    if ~isempty(ind)
        nanMap = (labelMap == ind);
    end


    blurMap = nanMap;
    if configS.showImgAnalFigs
        myImageSc(blurMap);
    end
else
    convMap = conv2(1.0*nanMap, templateMap, 'same');
    if configS.showImgAnalFigs
        myImageSc(convMap);
    end
    %myImageSc(convMap);
    % Find everythign that looks like a pit
    peaks = convMap > imdilate(convMap, [1 1 1; 1 0 1; 1 1 1]);
    peaks = 1.0*peaks .* convMap;
    if configS.showImgAnalFigs
        myImageSc(peaks);
    end
    %h = fspecial('gaussian', 8, 4);
    %myImageSc(h);
    %peaks = imfilter(peaks, h);
     if configS.showImgAnalFigs
        myImageSc(peaks);
    end
    [~,~,peakVals] = find(peaks);
    peakVals = sort(peakVals, 'descend');
    %disp(['   Max conv values:']);
    %peakVals(1:4)
    peakValCutoff = 0.8;
    if peakVals(1) < 0
        error('The peak is less than 0!');
    end
    peakVals(peakVals < peakValCutoff*max(peakVals)) = [];
    %peakVals
    x = zeros(1,length(peakVals));
    y = zeros(1,length(peakVals));
    for i=1:length(peakVals)
        %y(i)
        %x(i)
        %[valy valx] = find(peaks == peakVals(i), 1, 'first')
        [y(i),x(i)] = find(peaks == peakVals(i), 1, 'first');
        peaks(y(i), x(i)) = nan;
    end
    xc = round(length(nanMap(1,:))/2);
    yc = round(length(nanMap(:,1))/2);
    distToCenter = (x - xc).^2 + (y - yc).^2;
    [~,ind] = min(distToCenter(:));
    xBest = x(ind);
    yBest = y(ind);
    
    areaTemplateMap = sum(sum(1.0*(templateMap >= 1)));
    filledMask = imfill(1.0*(templateMap >= 1));
    areaFilledMask = sum(filledMask(:));
    strelPx = 1;
    while areaFilledMask < 1.2*areaTemplateMap
        filledMask = imclose(filledMask, strel('disk',strelPx));
        filledMask = imfill(filledMask);
        areaFilledMask = sum(filledMask(:));
        strelPx = strelPx + 1;
    end
    %shapeEdge = edge(filledMask);
    %shapeEdge = shapeEdge / sum(shapeEdge(:));
    SEout = strel('disk', round(configS.templateMaskOuter / scale));
    outMask = imdilate(filledMask, SEout) - filledMask;
    SEin = strel('square', round(configS.templateMaskInner / scale));
    inMask = xor(imerode(filledMask, SEin), filledMask);
    maskShape = outMask + inMask;
    if configS.showImgAnalFigs == 1
        myImageSc(filledMask);
        %myImageSc(shapeEdge);
        myImageSc(outMask);
        myImageSc(inMask);
        myImageSc(maskShape);
    end


    maskShape( ~any(maskShape,2), : ) = [];  %rows
    maskShape( :, ~any(maskShape,1) ) = [];  %columns

    if configS.showImgAnalFigs == 1
        myImageSc(maskShape);
    end

    mask = zeros(size(nanMap));
    mask(yBest,xBest) = 1;
    %myImageSc(mask);
    mask = conv2(mask, maskShape, 'same');
    %myImageSc(mask);
    %myImageSc(nanMap);
    blurMap = nanMap;
    %myImageSc(blurMap);
end

% An initial open operation gets rid of noise and removes connectivity of
% small shapes
%SE = strel('square', 3);
%openMap = imopen(nanMap, SE);
%if configS.showImgAnalFigs
%    myImageSc(openMap);
%end

%  - Gaussian blur (remove high frequency noise)
% BA - try smaller values
%blurFilter = fspecial('gaussian', 15*[1 1], 5);
%blurMap = imfilter(1.0*nanMap,blurFilter);
%blurMap = nanMap;
%if configS.showImgAnalFigs
%    myImageSc(blurMap);
%end

%  - Edge detection
edgeMap = edge(blurMap);
if ~isempty(mask)
    edgeMap = edgeMap .* mask;
end
if configS.showImgAnalFigs
    myImageSc(edgeMap);
end



end

