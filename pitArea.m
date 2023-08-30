function [ D1, boundaryWithoutOutliers ] = pitArea( boundary, configS, ...
    scale, doOutlierDetection, figHS )
%PITAREA Find the area of the pit
%   Detailed explanation goes here

ellipseBoundaryMaxDist = configS.ellipseOutlierMaxDist / scale;
ellispeOutlierN = configS.ellipseOutlierDiscardNPoints;

D1.a = nan;
D1.b = nan;
D1.phi = nan;
D1.X0 = nan;
D1.Y0 = nan;
D1.X0_in = nan;
D1.Y0_in = nan;
D1.long_axis = nan;
D1.short_axis = nan;
boundaryWithoutOutliers = [];

numPoints = [];
maxDist = [];
    
bx = boundary(:,2);
by = boundary(:,1);
if configS.showEllipseFigs == 1
    figure();
end

% Loop until there are only 5 points left; this is the minimum needed to
% define an ellipse. Discard all other points are being too erroneous (e.g.
% aggressive outlier detection). Could probably set a better condition
% considering residuals.
for i=1:length(bx)-10
    %  - Fit an ellipse to the data (least-squares method)
    E1 = fit_ellipse(bx,by);
    if ~isempty(E1.status);
        disp(['    **** ' E1.status ' ****']);
        D1.a = nan;
        D1.b = nan;
        D1.phi = nan;
        D1.X0 = nan;
        D1.Y0 = nan;
        D1.X0_in = nan;
        D1.Y0_in = nan;
        D1.long_axis = nan;
        D1.short_axis = nan;
        boundaryWithoutOutliers = [];
        return;
    end
    D1 = E1;
    % Sometimes a > b; othertimes a < b (a is always along the x-axis)
    D1.phi = -1*E1.phi;
    if E1.a < E1.b
        D1.a = E1.b;
        D1.b = E1.a;
        D1.phi = D1.phi + pi/2;
    end
    if D1.phi < 0
        D1.phi = D1.phi + pi;
    end
    if doOutlierDetection == 0
        boundaryWithoutOutliers = boundary;
        return;
    end
        

    %  - Outlier detection
    a = D1.a;
    b = D1.b;
    c = [D1.X0_in D1.Y0_in 0];
    u = [cos(D1.phi) sin(D1.phi) 0];
    v = [-sin(D1.phi) cos(D1.phi) 0];

    min_dist = distanceEllipsePoints([bx by zeros(size(bx))], a,b,c,u,v);

    [dists, distsI] = sort(min_dist);

    if configS.showEllipseOutlier == 1 || configS.showEllipseFigs == 1
        numPoints(i) = length(bx);
        maxDist(i) = dists(end) * scale;
    end
    
    % If all of the points are fairly close, as defined by
    % ellipseBoundaryMaxDist, then the ellipse is a pretty good fit, so
    % stop.
    if dists(end) < ellipseBoundaryMaxDist
        break;
    end
    
    % Remove the point with the largest residual
    bx(distsI(end-ellispeOutlierN+1:end)) = [];
    by(distsI(end-ellispeOutlierN+1:end)) = [];
    

end

if configS.showEllipseFigs == 1
    %[X, Y] = calcEllipse(E, 360);   
    %plot(X,Y, 'r-');
    %hold on;
    plot(boundary(:,2),boundary(:,1), 'rx');
    axis equal;
    hold on;
    plot(bx, by, 'bx');
    t = linspace(0,2*pi,360);
    x = D1.X0_in + D1.a*cos(t)*cos(D1.phi) - D1.b*sin(t)*sin(D1.phi);
    y = D1.Y0_in + D1.a*cos(t)*sin(D1.phi) + D1.b*sin(t)*cos(D1.phi);
    plot(x, y, 'k-');
    %pause(0.1);
end
if configS.showEllipseOutlier == 1
    set(0,'CurrentFigure',figHS.ellipseOutlier);
elseif configS.showEllipseFigs == 1
    figure();
end
if configS.showEllipseFigs == 1 || configS.showEllipseOutlier == 1
    plot(numPoints, maxDist, 'ko');
    xlabel('Number of boundary points');
    ylabel('Maximum distance from a boundary point to the ellipse');
    set(gca, 'XDir', 'Reverse');
end

boundaryWithoutOutliers = [bx by];

%area = pi*D1.a*D1.b;
end

