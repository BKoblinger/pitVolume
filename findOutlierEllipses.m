function [ outlierIndex, meanA, stdA, meanB, stdB, meanPhi ] = findOutlierEllipses( pitGeomS )
%FINDOUTLIERELLIPSES Summary of this function goes here
%   Detailed explanation goes here
phi = [pitGeomS.ellipsePhi];
a = [pitGeomS.ellipseA];
b = [pitGeomS.ellipseB];
    
outlierIndex = isnan(a);
prevNumOutliers = -1;
numOutliers = 0;
while numOutliers > prevNumOutliers
    prevNumOutliers = numOutliers;
    
    meanPhi = mean(phi(~isnan(a) & ~outlierIndex));
    %stdPhi = std([pitGeomS.ellipsePhi]);
    meanA = mean(a(~isnan(a) & ~outlierIndex));
    stdA = std(a(~isnan(a) & ~outlierIndex));
    meanB = mean(b(~isnan(a) & ~outlierIndex));
    stdB = std(b(~isnan(a) & ~outlierIndex));

    outlierIndex = (abs(a - meanA) > 2.0*stdA | abs(b - meanB) > 2.0*stdB) | isnan(a);
    numOutliers = sum(outlierIndex);
end

%outlier = (pitGeomS.ellipseA - meanA)^2 + (pitGeomS.ellipseB - meanB)^2;
end

