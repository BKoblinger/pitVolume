function [ avgElevation, sdFromPlane, stdg0, stdl0 ] = elevationPolishedSurface( xc, yc, ...
    configS, mapS, figHS )
%ELEVATIONPOLISHEDSURFACE Get the average elevation of the polished surface
%   Detailed explanation goes here

%  - Get appropriate radius boundaries from user

innerRadiusPolishedSurface = ...
    configS.polishedSurfaceInnerDiameter / mapS.scale / 2;
outerRadiusPolishedSurface = ...
    configS.polishedSurfaceOuterDiameter / mapS.scale / 2;

planeOutlierMaxDist = configS.planeOutlierMaxDist;
planeOutlierN = configS.planeOutlierDiscardNPoints;

circleMask = (mapS.y - yc).^2 + (mapS.x - xc).^2 <= (outerRadiusPolishedSurface)^2 & ...
             (mapS.y - yc).^2 + (mapS.x - xc).^2 > innerRadiusPolishedSurface^2;
circleMaskThreshold = (mapS.y - yc).^2 + (mapS.x - xc).^2 <= (outerRadiusPolishedSurface)^2 & ...
    (mapS.y - yc).^2 + (mapS.x - xc).^2 > innerRadiusPolishedSurface^2;

         
if configS.showPolishedElevFigs == 1
    myImageSc(mapS.mapZ.*circleMask);
end

[y, x] = find(circleMask);
z = mapS.mapZ(sub2ind(size(mapS.mapZ), y, x));
x = x(~isnan(z));
y = y(~isnan(z));
z = z(~isnan(z));


xorig = x;
yorig = y;
zorig = z;

if isempty(z)
    disp('"z" is empty. In "elevationPolishedSurface()"');
    avgElevation = nan;
    sdFromPlane = nan;
    stdg0 = nan;
    stdl0 = nan;
    return;
end


% If there might be grain edges
if configS.pitsNearGrainEdges == 1
    [yTh, xTh] = find(circleMaskThreshold);
    zTh = mapS.mapZ(sub2ind(size(mapS.mapZ), yTh, xTh));
    zTh = zTh(~isnan(z));
    % Try to threshold the polished mineral surface from the epoxy.
    histBinDZ = 1; % nm; must be an integer
    delZ = configS.epoxyGrainElevationDifference;
    gausHalfWidth = delZ/4;
    gausZ2 = -gausHalfWidth:histBinDZ:gausHalfWidth+delZ;
    %gausZ1 = -gausHalfWidth:histBinDZ:gausHalfWidth;
    c = gausHalfWidth*sqrt(2)/2;
    gaus2 = (gauss_distribution(1, 0, c, gausZ2) ...
        + gauss_distribution(1, delZ, c, gausZ2)) / ((c/sqrt(2))*sqrt(2*pi));
    %gaus1 = gauss_distribution(1, 0, c, gausZ1) / ((c/sqrt(2))*sqrt(2*pi));
    %figure();
    %plot(gausZ2, gaus2);
    %figure();
    %plot(gausZ1, gausZ1);
    
    
    
    [counts,centers] = hist(zTh, min(zTh):histBinDZ:max(zTh));
    filter2 = conv(counts, gaus2, 'same');
    %filter1 = conv(counts, gaus1, 'same');
    % Now smooth the data a bit
    windowSize = delZ/2;
    F = ones(1,windowSize)/windowSize;
    filter2smooth = conv(filter2,F, 'same');
    %filter1 = conv(filter1,F, 'same');
    
    
    %figure();
    %plot(centers, filter1);
    
    [peaks2,peakLocs2] = findpeaks(filter2smooth, 'MINPEAKHEIGHT', 10);
    %centers(peakLocs2)
    %[peaks1,peakLocs1] = findpeaks(filter1, 'MINPEAKHEIGHT', 10)
    
    if length(peaks2) >= 3
        disp('    Found epoxy and polished surface');
        [~,maxPeakInd] = max(peaks2);
        zCutoff = centers(peakLocs2(maxPeakInd));
        
        disp(['    Thresholding for grain versus epoxy at z = ' num2str(zCutoff)]);
        x(z < zCutoff) = [];
        y(z < zCutoff) = [];
        %zold = z;
        z(z < zCutoff) = [];
        
        epoxyMap = mapS.mapZ < zCutoff;
        
        
        if configS.showPolishedSurfaceHistogram == 1 ...
            || configS.showPolishedElevFigs == 1
            scrsz = get(0,'ScreenSize');
            figure('OuterPosition',[scrsz(3)/4 0 scrsz(3)/4 scrsz(4)/3], 'Visible', 'on');
            [countsThresh] = hist(z, centers);
            h = area(centers, counts, 'EdgeColor','none');
            set(h(1), 'FaceColor', [1 0 0]);
            hold on;
            h = area(centers, countsThresh, 'EdgeColor','none');
            set(h(1), 'FaceColor', [0 0 1]);
            plot(centers, filter2, 'k-', centers, filter2smooth, 'k-');
        end
        if configS.showPolishedElevFigs == 1    
            figure();
            myImageSc(epoxyMap);
            cm = flipud(brewermap(64,'RdBu'));
            colormap(cm);
            %scatter3(xorig, yorig, zorig, 2, 'r');
            %hold on;
            %scatter3(x, y, z, 2, 'k');
            
        end
        SE = strel('disk', 5);
        epoxyMap = imclose(epoxyMap, SE);
        if configS.showPolishedElevFigs == 1
            figure();
            myImageSc(epoxyMap);
            cm = flipud(brewermap(64,'RdBu'));
            colormap(cm);
        end
        grainMap = mapS.mapZ;
        grainMap(epoxyMap) = nan;
        if configS.showPolishedElevFigs == 1
            figure();
            myImageSc(grainMap);
            cm = flipud(brewermap(64,'RdBu'));
            colormap(cm);
        end
        [y, x] = find(circleMask);
        z = grainMap(sub2ind(size(epoxyMap), y, x));
        x = x(~isnan(z));
        y = y(~isnan(z));
        z = z(~isnan(z));
    else
      disp('   Did not find epoxy/polished surface division.');
    end
end
       
    
   
    %{
    f = fit(centers.',counts.','gauss2');

    if configS.showPolishedElevFigs == 1
        figure();
        hist(z, centers);
        hold on;
        plot(centers,gauss_distribution(f.a1, f.b1, f.c1, centers), 'r-');
        plot(centers,gauss_distribution(f.a2, f.b2, f.c2, centers), 'g-');
    end
    asd

    % Check to see if I did find a reasonably well defined epoxy peak. If so,
    % threshold, if not, don't.
    if f.a1 > 0 && f.a2 > 0 ...
            && (f.a1 > f.a2 && f.a1 < 10*f.a2) || (f.a2 > f.a1 && f.a2 < 10*f.a1)
        % Find the epoxy.
        bLo = f.b2;
        cLo = f.c2/sqrt(2);
        bHi = f.b1;
        cHi = f.c1/sqrt(2);
        if f.b1 < f.b2
            bLo = f.b1;
            cLo = f.c1/sqrt(2);
            bHi = f.b2;
            cHi = f.c2/sqrt(2);
        end   

        %Check to see if the two peaks overlap within two sigma. Only if there
        %are two very well resolved peaks, should data be eliminated
        loHighBound = bLo + 2*cLo;
        hiLowBound  = bHi - 2*cHi;
        if loHighBound < hiLowBound

            % check to see if there is a crossing point for the two gaussians
            % and only exclude up to the crossing.
            %zCross1 = (f.b2*f.c1^2 - f.b1*f.c2^2 + f.c1*f.c2*(f.b1^2 - 2*f.b1*f.b2 + f.b2^2 - log(f.a1/f.a2)*f.c1^2 + log(f.a1/f.a2)*f.c2^2)^(1/2))/(f.c1^2 - f.c2^2)
            %zCross2 = -(f.b1*f.c2^2 - f.b2*f.c1^2 + f.c1*f.c2*(f.b1^2 - 2*f.b1*f.b2 + f.b2^2 - log(f.a1/f.a2)*f.c1^2 + log(f.a1/f.a2)*f.c2^2)^(1/2))/(f.c1^2 - f.c2^2)
            %if gauss_distribution(f.a1, f.b1, f.c1, zCross1) > 1 && gauss_distribution(f.a2, f.b2, f.c2, zCross1)
            %    zMax = zCross1
            %elseif gauss_distribution(f.a1, f.b1, f.c1, zCross2) > 1 && gauss_distribution(f.a2, f.b2, f.c2, zCross2)
            %    zMax = zCross2
            %else
                %zMax = (z0 + 3.0*c)
            %end
            disp('    Thresholding for grain versus epoxy.');
            x(z < loHighBound) = [];
            y(z < loHighBound) = [];
            z(z < loHighBound) = [];
        end
    end

    if configS.showPolishedElevFigs == 1
        figure();
        hist(z, centers);
        hold on;
        plot(centers,gauss_distribution(f.a1, f.b1, f.c1, centers), 'r-');
        plot(centers,gauss_distribution(f.a2, f.b2, f.c2, centers), 'g-');
    end
end
    %}


%a = [];
%b = [];
%c = [];
%d = [];
sdFromPlane = [];
%figure();
for i=1:length(z)-3
    % Fit the data to a plane
    % z(x,y) = p00 + p10*x + p01*y;
    planefitobject = fit([x,y],z,'poly11');

    coef = coeffvalues(planefitobject);
    %confint(planefitobject); % confidence intervals

    %  - Outlier detection
    a = coef(2);
    b = coef(3);
    c = -1;
    d = coef(1);
    D = (a*x + b*y + c*z + d) ./ sqrt(a^2 + b^2 + c^2);
    
    [dist, distI] = sort(abs(D));
    if configS.showPlaneOutlier == 1
        maxDist(i) = dist(end);
        numPoints(i) = length(z);
    end
    
    %plot(numPoints(i), maxDist(i), 'kx');
    %hold on;
    %pause(0.001);
    

    % Get rid of the N points farthest away from the plane.
    if length(distI) > 2*planeOutlierN
        x(distI(end-planeOutlierN+1:end)) = [];
        y(distI(end-planeOutlierN+1:end)) = [];
        z(distI(end-planeOutlierN+1:end)) = [];
    else
        disp(['    **** There are only ' num2str(length(distI)) ' points left; will not discard them; elevationPolishedSurface.m']); 
        break;
    end
    
    if dist(end) < planeOutlierMaxDist
        break;
    end
end

D = (a*x + b*y + c*z + d) ./ sqrt(a^2 + b^2 + c^2);
sdReducedZ = stdZeroMean(D);
avgElevation = coef(1) + coef(2)*xc + coef(3)*yc;
D = (a*xorig + b*yorig + c*zorig + d) ./ sqrt(a^2 + b^2 + c^2);
sdFromPlane1 = stdZeroMean(D);
% Now get rid of outliers greater than 2 sigma
isTooFar = abs(D) > 2*sdFromPlane1;
%x2 = xorig;
%y2 = yorig;
%z2 = zorig;
%x2(isTooFar) = [];
%y2(isTooFar) = [];
%z2(isTooFar) = [];
sdFromPlane = stdZeroMean(D(~isTooFar));

Dg0 = D(D>=0);
Dl0 = D(D<0);
if length(Dg0) + length(Dl0) ~= length(D)
    error('Well this is weird.');
end
stdg0 = stdZeroMean(Dg0);
stdl0 = stdZeroMean(Dl0);
Dg0 = Dg0(abs(Dg0) <= 2*stdg0);
Dl0 = Dl0(abs(Dl0) <= 2*stdl0);
stdg0 = stdZeroMean(Dg0);
stdl0 = stdZeroMean(Dl0);




if configS.showPolishedElevFigs == 1
    figure();
    scatter3(xorig, yorig, zorig, 2, 'r');
    hold on;
    scatter3(x, y, z, 2, 'k');
    hold on;
    grid on;

    [yp, xp] = size(mapS.mapZ);
    [XX, YY] = meshgrid([1 xp], [1 yp]);
    Z = coef(1) + coef(2)*XX + coef(3)*YY;
    reorder = [1 2  4 3];
    patch(XX(reorder), YY(reorder), Z(reorder), 'b');
    alpha(0.3);
    view(40,35)


    xlabel('x');
    ylabel('y');
    zlabel('z');
end

if configS.showPlaneOutlier == 1 
    set(0, 'CurrentFigure', figHS.planeOutlier);
elseif configS.showPolishedElevFigs == 1
    figure();
end
if configS.showPlaneOutlier == 1 || configS.showPolishedElevFigs == 1
    plot(numPoints, maxDist, 'ko');
    xlabel('Number of points in plane-fitting calculation');
    ylabel('Distance from the plane to the farthest point in calculation (nm)');
    set(gca, 'XDir', 'Reverse');
end

    
if ~isempty(figHS.summary)
    set(0, 'CurrentFigure', figHS.summary)
    subplot(3,2,3,'replace');

    if ~isempty(z)
        sdz = std(z);
        medianz = median(z);
        pitImageSc(mapS, medianz + 2*sdz*[-1 1]);
    else
        pitImageSc(mapS);
    end
    
    hold on;
    t = linspace(0,2*pi,360);
    % Plot the "polished surface" circles.
    xCircIn  = mapS.scale * (xc + innerRadiusPolishedSurface*cos(t));
    yCircIn  = mapS.scale * (yc + innerRadiusPolishedSurface*sin(t));
    xCircOut = mapS.scale * (xc + outerRadiusPolishedSurface*cos(t));
    yCircOut = mapS.scale * (yc + outerRadiusPolishedSurface*sin(t));

    plot(xCircIn, yCircIn, 'k-');
    plot(xCircOut, yCircOut, 'k-');
    title('Polished surface (nm)');


    subplot(3,2,4,'replace');
    
    plane = coef(1) + coef(2)*mapS.x + coef(3)*mapS.y;
    dMapS = mapS;
    dMapS.mapZ = dMapS.mapZ - plane;
    pitImageSc(dMapS, 2*sdFromPlane*[-1 1]);
    
    %hold on;
    %contour(dImg, [0 0], 'b-');
    %contour(dImg, [-1*sd sd], 'c-');
    
    hold on;
    plot(xCircIn, yCircIn, 'k-');
    plot(xCircOut, yCircOut, 'k-');
    title('Polished surface corrected for tilt (\Delta nm)');
end

end


function f = gauss_distribution(a1, b1, c1, x)
f = a1*exp(-((x-b1)/c1).^2);
end

function stdev = stdZeroMean(x)
stdev = sqrt(sum(x.^2)/(length(x)-1));
end



