function [ boundaryOut ] = reduceBoundary( mapS, boundary, xc, yc, keepPoints )
%REDUCEBOUNDARY Summary of this function goes here
%   Detailed explanation goes here

%disp('WARNING: the parameters keepPoints has not been fully tested. It is best to set keepPoints = all');

lims = size(mapS.mapZ);
xlo = 1;
xhi = lims(2);
ylo = 1;
yhi = lims(1);

boundaryOut = zeros(size(boundary));

boundaryM = zeros(size(mapS.mapZ));
boundaryM(sub2ind(size(boundaryM),boundary(:,1),boundary(:,2))) = 1;
%myImageSc(boundaryM);
for i=1:length(boundary(:,1))
    yp = boundary(i,1);
    xp = boundary(i,2);
    %{
    m = (yp - yc)/(xp - xc);
    b = yp - m*xp;
    
    if xp > xc
        x2 = xhi;
    else
        x2 = xlo;
    end
    
    y2 = m*x2 + b;
    if y2 > yhi || y2 < ylo
        if yp > yc
            y2 = yhi;
        else
            y2 = ylo;
        end
        x2 = (y2 - b)/m;
    end
    %}
    
    %[x1, y1; x2, y2]
    %[cx,cy,c] = improfile(boundaryM,[xc x2] ,[yc y2]);
    [c, cy, cx] = bresenham(boundaryM,[yc xc; yp xp], 0);

    %[myline,mycoords,outmat,X,Y] = bresenham(mymat,mycoordinates,dispFlag)
    %myImageSc(boundaryM);
    %hold on;
    %plot([xc x2] ,[yc y2]);
    %hold on;
    %plot(xc, yc, 'gx');
    %plot(xp, yp, 'rx');
    
    %figure();
    %plot(c);
    
    
    
    cx = cx(c ~= 0);
    cy = cy(c ~= 0);
    
    dist = (cx - xc).^2 + (cy - yc).^2;
    if keepPoints == 1
        [~,minInd] = min(dist);
    elseif keepPoints == 2
        [~,minInd] = max(dist);
    else
        error('The variable keepPoints was not set correctly.');
    end
    distBoundary = (cx(minInd) - boundary(:,2)).^2 + (cy(minInd) - boundary(:,1)).^2;
    [~,minInd] = min(distBoundary);
    boundaryOut(i, :) = boundary(minInd, :);
end

boundaryOut = unique(boundaryOut,'rows');
%boundaryM = zeros(size(mapS.mapZ));
%boundaryM(sub2ind(size(boundaryM),boundaryOut(:,1),boundaryOut(:,2))) = 1;
%myImageSc(boundaryM);


end

