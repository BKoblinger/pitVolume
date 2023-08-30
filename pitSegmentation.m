function [ boundary ] = pitSegmentation( mapS, convMap, configS, figHS )
%PITSEGMENTATION Summary of this function goes here
%   Detailed explanation goes here

approxPitRadiusPixels = configS.pitDiameter / mapS.scale / 2;
approxPitFloorRadiusPixels = configS.pitFloorDiameter / mapS.scale / 2;
houghDMaskInner = configS.houghMaskInner / mapS.scale;
houghDMaskOuter = configS.houghMaskOuter/ mapS.scale;


nanMap = isnan(mapS.mapZ);
if configS.showImgAnalFigs == 1
    myImageSc(nanMap);
end


[edgeMap, templateMask] = deNoise(nanMap, mapS.scale, convMap, approxPitRadiusPixels, ...
    approxPitFloorRadiusPixels, configS);
    
if isempty(convMap)
    % Fit a Hough circle. Use this as the approximate location of the ellipse.
    [centers, radii] = imfindcircles(edgeMap,...
        round([0.9*approxPitRadiusPixels 1.1*approxPitRadiusPixels]),...
        'Sensitivity',0.99);

    if configS.showImgAnalFigs == 1
        if ~isempty(radii)
            viscircles(centers(1,:), radii(1),'EdgeColor','b');
        end
        if length(radii) >= 2
            viscircles(centers(2,:), radii(2),'EdgeColor','c');
        end
        %viscircles(centers(3,:), radii(3),'EdgeColor','r');
    end

    %dRadius = 0.15*radii(1);

    % In the previous step we fit a circle to the data and now we turn the 
    % cirlce into a disc and make a mask to get rid of data "far away" from the
    % ellipse. If the pits are elliptical, we need to use a larger mask, but
    % this increases the effect of non-cicular "defects."
    dInRadius = houghDMaskInner;
    dOutRadius = houghDMaskOuter;
    
    if ~isempty(centers)

        houghMask = (mapS.y - centers(1,2)).^2 + (mapS.x - centers(1,1)).^2 <= (radii(1) + dOutRadius)^2 & ...
                     (mapS.y - centers(1,2)).^2 + (mapS.x - centers(1,1)).^2 >  (radii(1) - dInRadius)^2;

        
        if configS.showHoughMask == 1
            set(0, 'CurrentFigure', figHS.houghMask);
        elseif configS.showImgAnalFigs == 1 
            figure();
        end
        if configS.showHoughMask == 1 || configS.showImgAnalFigs
            houghMapS = mapS;
            houghMapS.mapZ = edgeMap + houghMask;
            pitImageSc(houghMapS);
            hold on;
            t = linspace(0,2*pi,360);
            % Plot the "polished surface" circles.
            xCirc  = centers(1,1) + radii(1)*cos(t);
            yCirc  = centers(1,2) + radii(1)*sin(t);
            plot(xCirc*mapS.scale, yCirc*mapS.scale, 'g-');
            title('Hough mask extent');
            pause(0.1);
        end

        edgeMap = edgeMap .* houghMask;
    else
        houghMask = [];
        disp('    Could not fit a circle to the edges. Something has probably gone wrong in pitSegmentation.');
    end
end
    
if configS.findRim == 1
    if isempty(convMap)
        [y, x] = find(edgeMap);
        boundary1 = [y x];

        ellipseS = pitArea(boundary1, configS, mapS.scale, 0);
        a = ellipseS.a;
        b = ellipseS.b;
        phi = ellipseS.phi;
        xc = ellipseS.X0_in;
        yc = ellipseS.Y0_in;
        thresholdShape = ...
            (((mapS.x-(xc))*cos(phi)+(mapS.y-(yc))*sin(phi)).^2/(a)^2) + ...
            (((mapS.x-(xc))*sin(phi)-(mapS.y-(yc))*cos(phi)).^2/(b)^2) <= 1;  
        SE = strel('disk', round(configS.houghMaskOuter/mapS.scale));
        thresholdMask = imdilate(thresholdShape, SE);
        thresholdMask = thresholdMask - thresholdShape;
    else
        thresholdMask = templateMask;
    end
    

    thresholdVals = mapS.mapZ .* thresholdMask;
    if configS.showImgAnalFigs == 1
        myImageSc(thresholdVals);
    end
    thresholdVals = thresholdVals(~isnan(thresholdVals));
    thresholdVals = thresholdVals(thresholdVals ~= 0);

    sdThreshold = std(thresholdVals);
    m0Threshold = median(thresholdVals);

    rim = mapS.mapZ;
    rim(rim < m0Threshold + configS.findRimSigma*sdThreshold)  = nan;
    if configS.showImgAnalFigs == 1
        myImageSc(rim);
    end

    if configS.showImgAnalFigs == 1
        histBinDZ = 1; % nm; must be an integer
        [~,centers] = hist(thresholdVals, min(thresholdVals):histBinDZ:max(thresholdVals));
        %windowSize = 100;
        %F = ones(1,windowSize)/windowSize;
        %smoothed = conv(counts,F, 'same');
        figure();
        hist(thresholdVals, centers);
        hold on;
        hist(rim(:), centers);
    end

    rim(isnan(rim)) = 0;

    edgeMap = edge(rim);
    if isempty(convMap)
     edgeMap = edgeMap .* houghMask;
    else
        edgeMap = edgeMap .* templateMask;
    end
    if configS.showImgAnalFigs == 1
        myImageSc(edgeMap);
    end

    %edgeMap = deNoise(rim, approxPitRadiusPixels, ...
    %    approxPitFloorRadiusPixels, configS.showImgAnalFigs);


    %[centers, radii] = imfindcircles(edgeMap,...
    %    round([0.9*approxPitRadiusPixels 1.1*approxPitRadiusPixels]),...
    %    'Sensitivity',0.99);

    %if configS.showImgAnalFigs == 1
    %    if ~isempty(centers)
    %        viscircles(centers(1,:), radii(1),'EdgeColor','b');
    %        %viscircles(centers(2,:), radii(2),'EdgeColor','c');
    %        %viscircles(centers(3,:), radii(3),'EdgeColor','r');
    %    end
    %end
end
        

if configS.showImgAnalFigs == 1
    myImageSc(edgeMap);
end


[y, x] = find(edgeMap);
if ~isempty(y)
    boundary = [y x];
else
    boundary = [];
end


end

%{
if 0
    SE = strel('disk', 5);
    edgeMap = imclose(edgeMap, SE);
    if configS.showImgAnalFigs == 1
        myImageSc(edgeMap);
    end

    %C = imclearborder(edgeMap);
    C = edgeMap;

    % Try to make the outline a solid shape. This gets rid of small gaps in the
    % ellipse outline. Big holes will not be filled.

    SE = strel('square', 3);
    %C = imclose(C, SE);
    C = imdilate(C, SE);
    if configS.showImgAnalFigs == 1
        myImageSc(C);
    end
    %C = imfill(C, 'holes');
    %if configS.showImgAnalFigs == 1
    %    myImageSc(C);
    %end
    C = imerode(C, SE);
    if configS.showImgAnalFigs == 1
        myImageSc(C);
    end


    %C = imfill(C,'holes');
    %if configS.showImgAnalFigs == 1
    %    myImageSc(edgeMap);
    %end


    %[B,L,N,A] = bwboundaries(C,8,'noholes');
    %[B,L,~,A] = bwboundaries(C,8,'noholes');
    [B,L,~,A] = bwboundaries(C,8);
    if configS.showImgAnalFigs == 1
        colors=['b' 'g' 'r' 'c' 'm'];
        myImageSc(mapS.mapZ); hold on;
        for k=1:length(B)
            boundary = B{k};
            %boundarySize(k) = (max(boundary(:,1)) - min(boundary(:,1))) + ...
            %                  (max(boundary(:,2)) - min(boundary(:,2)));
            cidx = mod(k,length(colors))+1;
            plot(boundary(:,2), boundary(:,1),...
                 colors(cidx),'LineWidth',2);
            %randomize text position for better visibility
            rndRow = ceil(length(boundary)/(mod(rand*k,7)+1));
            col = boundary(rndRow,2); row = boundary(rndRow,1);
            h = text(col+1, row-1, num2str(L(row,col)));
            set(h,'Color',colors(cidx),...
                'FontSize',14,'FontWeight','bold');
        end
    end

    % enclosing_boundary  = find(A(m,:));
    % enclosed_boundaries = find(A(:,m));


    % Include all boundaries inside the circle mask (from the Hough transform)


    lenBoundary = zeros(1,length(B));
    for i=1:length(B)  

        %if isempty(circleMask)
            lenBoundary(i) = length(B{i});
        %elseif ~any(A(i,:)) ...
        %        && sum(circleMask(sub2ind(size(circleMask), B{i}(:,1), B{i}(:,2)))) > 10
            % if ~any(A(i,:)) % It is enclosed by nothing
            % if sum(circleMask(sub2ind(size(circleMask), boundary(1,:), boundary(2,:)))) > 10
            %    There are several pixels 
        %    lenBoundary(i) = length(B{i});
        %else
        %    % It is enclosed by another boundary. Discard.
        %    lenBoundary(i) = 0;
        %end
        %enclosing_boundary  = find(A(i,:));

    end


    if configS.findRim == 0
        [~,ind] = find(lenBoundary > 0.1*pi*(approxPitRadiusPixels));
    else
        [~,ind] = find(lenBoundary > 0);
    end
    BB = B(ind);

    %%{
    %insiders = zeros(1,length(BB));
    %% My code to get rid of boundaries more or less inside other boundaries
    %for i=1:length(BB)
    %    bA = BB{i};
    %    for j=1:length(BB)
    %        bB = BB{j};
    %        if max(bA(:,1)) > max(bB(:,1)) && ...
    %                min(bA(:,1)) < min(bB(:,1)) && ...
    %                max(bA(:,2)) > max(bB(:,2)) && ...
    %                min(bA(:,2)) < min(bB(:,2))
    %            insiders(j) = 1;       
    %        end    
    %    end
    %end
    %BB = BB(~insiders);
    %%%


    if configS.showImgAnalFigs == 1
        colors=['b' 'g' 'r' 'c' 'm'];
        myImageSc(mapS.mapZ); hold on;
        for k=1:length(BB)
            boundary = BB{k};
            %boundarySize(k) = (max(boundary(:,1)) - min(boundary(:,1))) + ...
            %                  (max(boundary(:,2)) - min(boundary(:,2)));
            cidx = mod(k,length(colors))+1;
            plot(boundary(:,2), boundary(:,1),...
                 colors(cidx),'LineWidth',2);
        end
    end


    boundary = [];
    for i=1:length(BB)
        boundary = [boundary; BB{i}];
    end
    end
%}


%{
% Threshold to ignore the condensate rim.
    % Only do this for polished grains
    if configS.isPolished == 1 

        %[colsInImage, rowsInImage] = meshgrid(1:length(Map(1,:)), 1:length(Map(:,1)));
        thresholdMask = (mapS.y - centers(1,2)).^2 + (mapS.x - centers(1,1)).^2 <= (radii(1) + 2*dOutRadius)^2 & ...
                       (mapS.y - centers(1,2)).^2 + (mapS.x - centers(1,1)).^2 >  (radii(1))^2;

        thresholdVals = mapS.mapZ .* thresholdMask;
        thresholdVals = thresholdVals(~isnan(thresholdVals));
        thresholdVals = thresholdVals(thresholdVals ~= 0);

        sdThreshold = std(thresholdVals);
        m0Threshold = median(thresholdVals);

        if configS.showImgAnalFigs == 1
            figure();
            hist(thresholdVals, 300);
        end

        thresholdBoolMap = mapS.mapZ < (m0Threshold - 2.0*sdThreshold);
        thresholdMap = mapS.mapZ;
        thresholdMap(thresholdBoolMap == 1) = nan;
        %thresholdMap(thresholdBoolMap == 0) = Map(thresholdBoolMap == 0);

        nanMap = isnan(thresholdMap);

        if configS.showImgAnalFigs == 1
            myImageSc(thresholdMap);
            myImageSc(nanMap);
        end

        edgeMap = deNoise(nanMap, approxPitRadiusPixels, ...
            approxPitFloorRadiusPixels, configS.showImgAnalFigs);


        [centers, radii] = imfindcircles(edgeMap,...
            round([0.9*approxPitRadiusPixels 1.1*approxPitRadiusPixels]),...
            'Sensitivity',0.99);

        if configS.showImgAnalFigs == 1
            viscircles(centers(1,:), radii(1),'EdgeColor','b');
            viscircles(centers(2,:), radii(2),'EdgeColor','c');
            %viscircles(centers(3,:), radii(3),'EdgeColor','r');
        end
%}

