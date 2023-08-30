function plotPitEllipse(mapS, boundary, pitGeomS, scale, figHS)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

set(0,'CurrentFigure', figHS.summary);
%subplotH = 
subplot(3,2,1,'replace');

if ~isempty(boundary)
    borderMap = zeros(size(mapS.mapZ));
    borderMap(sub2ind(size(mapS.mapZ), boundary(:,1), boundary(:,2))) = 1;
    SE = strel('square', 15);
    borderMap = imdilate(borderMap, SE);
    borderMap(borderMap == 0) = nan;
    m0 = mapS.mapZ .* borderMap;
    m0 = m0(~isnan(m0));
    mm0 = median(m0);
    sdm0 = std(m0);

    pitImageSc(mapS, mm0 + 2*sdm0*[-1 1]);
else
    pitImageSc(mapS);
end

if ~isnan(pitGeomS.ellipseA)

    hold on;

    t = linspace(0,2*pi,360);
    xE = (scale*pitGeomS.ellipseXc + pitGeomS.ellipseA*cos(t)*cos(pitGeomS.ellipsePhi) - pitGeomS.ellipseB*sin(t)*sin(pitGeomS.ellipsePhi));
    yE = (scale*pitGeomS.ellipseYc + pitGeomS.ellipseA*cos(t)*sin(pitGeomS.ellipsePhi) + pitGeomS.ellipseB*sin(t)*cos(pitGeomS.ellipsePhi));
    plot(xE, yE, 'g-');

    plot([scale*pitGeomS.ellipseXc-pitGeomS.ellipseA*cos(pitGeomS.ellipsePhi) ...
        scale*pitGeomS.ellipseXc+pitGeomS.ellipseA*cos(pitGeomS.ellipsePhi)], ...
        [scale*pitGeomS.ellipseYc-pitGeomS.ellipseA*sin(pitGeomS.ellipsePhi) ...
        scale*pitGeomS.ellipseYc+pitGeomS.ellipseA*sin(pitGeomS.ellipsePhi)], 'b-');
    plot([scale*pitGeomS.ellipseXc+pitGeomS.ellipseB*sin(pitGeomS.ellipsePhi)...
        scale*pitGeomS.ellipseXc-pitGeomS.ellipseB*sin(pitGeomS.ellipsePhi)], ...
        [scale*pitGeomS.ellipseYc-pitGeomS.ellipseB*cos(pitGeomS.ellipsePhi) ...
        scale*pitGeomS.ellipseYc+pitGeomS.ellipseB*cos(pitGeomS.ellipsePhi)], 'g-');

    % Plot the "polished surface" circles.
    %xCircIn  = D.X0_in + innerRadius*cos(t);
    %yCircIn  = D.Y0_in + innerRadius*sin(t);
    %xCircOut = D.X0_in + outerRadius*cos(t);
    %yCircOut = D.Y0_in + outerRadius*sin(t);

    %plot(h, xCircIn, yCircIn, 'k-');
    %plot(h, xCircOut, yCircOut, 'k-');

    %dx = 0.1*max(xE);
    dy = 0.1*max(yE);
    yylim = [min(yE)-dy max(yE)+dy];
    xxlim = mean([min(xE) max(xE)]) + (4/3)*(yylim(2)-yylim(1)) * [-0.5 0.5];

    if isnan(yylim(1)) || isnan(yylim(2))
        yylim(1) = min(mapS.y(:)) * mapS.scale;
        yylim(2) = max(mapS.y(:)) * mapS.scale;
    end
    if isnan(xxlim(1)) || isnan(xxlim(2))
        xxlim(1) = min(mapS.x(:))* mapS.scale;
        xxlim(2) = max(mapS.x(:))* mapS.scale;
    end

    ylim(yylim);
    xlim(xxlim);
end

title('Pit edge and fitted ellipse');

end

