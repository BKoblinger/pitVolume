function [ medianElev, sdElev ] = pitFloorElevationSimple( xc, yc, mapS, configS, figHS )
%PITELEVATION Summary of this function goes here
%   Detailed explanation goes here

%  ** Get the average elevation in the centre of the pit
%  - User input (in config file) for diameter
%  - Get average elevation in circle
%pitFloorMeasRadius = configS.pitFloorMeasuredDiameter / mapS.scale / 2;
pitFloorMeasRadius = 18 / mapS.scale / 2;
 
%[colsInImage, rowsInImage] = meshgrid(1:length(Map(1,:)), 1:length(Map(:,1)));
circleMask = (mapS.y - yc).^2 + (mapS.x - xc).^2 < (pitFloorMeasRadius)^2;

elev = mapS.mapZ .* circleMask;

elev = elev(circleMask ~= 0);
elev = elev(~isnan(elev));
medianElev = median(elev);
sdElev = std(elev);

if isempty(elev)
    medianElev = nan;
    sdElev = nan;
end

if ~isempty(figHS.summary)
    set(0, 'CurrentFigure', figHS.summary);
    subplot(3,2,5,'replace');

    maxElev = max(max(elev(elev ~= 0)));
    minElev = min(min(elev(elev ~= 0)));
    dElev = (maxElev-minElev);
    if ~isempty(elev)
        pitImageSc(mapS, [minElev-0.5*dElev maxElev+0.1*dElev]);
    else
        pitImageSc(mapS);
    end
    if ~isnan(xc)
        xlim([xc - 4*pitFloorMeasRadius xc + 4*pitFloorMeasRadius] * mapS.scale);
        ylim([yc - 3*pitFloorMeasRadius yc + 3*pitFloorMeasRadius] * mapS.scale);
    end
    
    
    hold on;
    t = linspace(0,2*pi,360);
    x = xc + pitFloorMeasRadius*cos(t);
    y = yc + pitFloorMeasRadius*sin(t);
    plot(mapS.scale * x, mapS.scale * y, 'c-');
    title('Pit floor with analysis circle (nm)');
    
   
    h = subplot(3,2,6,'replace');
    dElev = 0.01E4;

    if ~isempty(elev)
        hist(elev, min(elev)-dElev/2:dElev:max(elev)+dElev/2);%, 'FaceColor', [1 0 0])
        %histogram(elev); , 'FaceColor', [1 0 0]
        hold on;
        vline(h, medianElev, 'r-');
        title('Floor elevations within analysis circle');
    end
end

end

