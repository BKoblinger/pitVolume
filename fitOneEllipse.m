function [ pitGeomS, boundary, borders, boundaryWithoutOutliers, isFit ] ...
    = fitOneEllipse( mapS, convMap, label, fileName, configS, showSumFigs, figHS, figHSprev  )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

pitGeomS = nullPitGeomS();
pitGeomS.file = fileName;
pitGeomS.label = label;

fprintf(1, ['Beginning ellipse fitting on file: "' pitGeomS.label '": ' pitGeomS.file]);

if showSumFigs
    plotOverviewMap( pitGeomS.label, figHS, mapS, 1 );
end

%% Find the ellipse border/boundar/outline {(x,y}) then fit an ellipse to it
fprintf(1, ' ... ');
startTime = tic();
%  ** Find the area of the pit
%  - Segmentation to get the edge of the pit as a matrix of (x,y) points
[ pitGeomS, boundary, borders, boundaryWithoutOutliers, isFit  ] ...
    = firstEllipseFitting( mapS, convMap, pitGeomS, configS, figHS );
ellapsedTime = toc(startTime);
fprintf(1, '%.2f s\n', ellapsedTime);

if showSumFigs
    plotEllipseFit(borders, boundary, ...
        boundaryWithoutOutliers, pitGeomS, mapS, figHS);
    plotPitEllipse(mapS, boundary, pitGeomS, mapS.scale, figHS);
    displayFigures( showSumFigs, pitGeomS.label, 1, figHS, figHSprev );
end

end

