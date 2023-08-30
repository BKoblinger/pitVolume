function [  ] = plotEllipseFit( border, boundary, boundaryWithoutOutliers, pitGeomS, mapS, figHS )
%PLOTELLIPSEFIT Summary of this function goes here
%   Detailed explanation goes here

set(0, 'CurrentFigure', figHS.summary);
%subplotH = 
subplot(3,2,2,'replace');
%subplot(subplotH);

if ~isempty(border)
    plot(mapS.scale*border(:,2),mapS.scale*border(:,1), 'gx', 'MarkerSize', 1);
else
    plot(-100,-100);
end
axis equal;
hold on;
if ~isempty(boundary)
    plot(mapS.scale*boundary(:,2),mapS.scale*boundary(:,1), 'rx', 'MarkerSize', 1);
end

if ~isempty(boundaryWithoutOutliers)
    plot(mapS.scale*boundaryWithoutOutliers(:,1), mapS.scale*boundaryWithoutOutliers(:,2), 'ko', 'MarkerSize', 1);
end
t = linspace(0,2*pi,360);
xE = mapS.scale*pitGeomS.ellipseXc + pitGeomS.ellipseA*cos(t)*cos(pitGeomS.ellipsePhi) - pitGeomS.ellipseB*sin(t)*sin(pitGeomS.ellipsePhi);
yE = mapS.scale*pitGeomS.ellipseYc + pitGeomS.ellipseA*cos(t)*sin(pitGeomS.ellipsePhi) + pitGeomS.ellipseB*sin(t)*cos(pitGeomS.ellipsePhi);

dy = 0.1*max(yE);
yylim = [min(yE)-dy max(yE)+dy];
xxlim = mean([min(xE) max(xE)]) + (4/3)*(yylim(2)-yylim(1)) * [-0.5 0.5];

if isnan(xxlim) + isnan(yylim) == 0
    ylim(yylim);
    xlim(xxlim);
end

%x = x * configS.mapS.scale;
%y = y * configS.mapS.scale;
plot(xE, yE, 'k-');
title('Pit edge and fitted ellipse');


end

