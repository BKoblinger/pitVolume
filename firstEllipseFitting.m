function [ pitGeomS, boundary, borders, boundaryWithoutOutliers, isFit  ] ...
    = firstEllipseFitting( mapS, convMap, pitGeomS, configS, figHS )
%FIRSTELLIPSEFITTING Summary of this function goes here
%   Detailed explanation goes here

%pitGeomS = nullPitGeomS();

boundary = [];
boundaryWithoutOutliers = [];
borders = pitSegmentation(mapS, convMap, configS, figHS);
        
isFit = isEllipseFit(borders, 'Could not find the ellipse border.');
if ~isFit
    return;
end

%if configS.isPolished == 1
%    configS.keepPoints = 'outside';
%else
%    keepPoints = 'inside';
%end

if configS.keepPoints >= 1
    ellipseS = pitArea(borders, configS, mapS.scale, 0);
    isFit = isEllipseFit(ellipseS, 'Could not fit an ellipse to the border 1.');
    if ~isFit
        return;
    end        

    boundary = reduceBoundary(mapS, borders, ellipseS.X0_in, ellipseS.Y0_in, configS.keepPoints);
else
    boundary = borders;
end

% Find the ellipse area in pixels and the ellipse parameters
[ellipseS, boundaryWithoutOutliers] = pitArea(boundary, configS, mapS.scale, 1, figHS);
isFit = isEllipseFit(ellipseS, 'Could not fit an ellipse to the border 2.');
if ~isFit
    return;
end

pitGeomS.ellipseA = ellipseS.a * (mapS.scale); % in um
pitGeomS.ellipseB = ellipseS.b * (mapS.scale); % in um
pitGeomS.ellipsePhi = ellipseS.phi; % in rad
pitGeomS.ellipseXc = ellipseS.X0_in; % in pixels!
pitGeomS.ellipseYc = ellipseS.Y0_in; % in pixels!
pitGeomS.ellipseArea = pi * pitGeomS.ellipseA * pitGeomS.ellipseB;  % in um^2

end

