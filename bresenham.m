function [myline,X,Y] = bresenham(mymat,mycoordinates,dispFlag)

% BRESENHAM: Generate a line profile of a 2d image 
%            using Bresenham's algorithm
% [myline,mycoords] = bresenham(mymat,mycoordinates,dispFlag)
%
% - For a demo purpose, try >> bresenham();
%
% - mymat is an input image matrix.
%
% - mycoordinates is coordinate of the form: [x1, y1; x2, y2]
%   which can be obtained from ginput function
%
% - dispFlag will show the image with a line if it is 1
%
% - myline is the output line
%
% - mycoords is the same as mycoordinates if provided. 
%            if not it will be the output from ginput() 
% Author: N. Chattrapiban
%
% Ref: nprotech: Chackrit Sangkaew; Citec
% Ref: http://en.wikipedia.org/wiki/Bresenham's_line_algorithm
% 
% See also: tut_line_algorithm

if nargin < 1, % for demo purpose
    pxl = 20;
    mymat = 1:pxl^2;
    mymat = reshape(mymat,pxl,pxl);
    disp('This is a demo.')
end

if nargin < 2, % if no coordinate provided
    imagesc(mymat); axis image;
    disp('Click two points on the image.')
    %[mycoordinates(1:2),mycoordinates(3:4)] = ginput(2);
    mycoordinates = ginput(2);
end

if nargin < 3, dispFlag = 1; end

%BMH
%outmat = mymat;
mycoords = mycoordinates;

x = round(mycoords(:,1));
y = round(mycoords(:,2));
steep = (abs(y(2)-y(1)) > abs(x(2)-x(1)));

if steep, [x,y] = swap(x,y); end

if x(1)>x(2), 
    [x(1),x(2)] = swap(x(1),x(2));
    [y(1),y(2)] = swap(y(1),y(2));
end

delx = x(2)-x(1);
dely = abs(y(2)-y(1));
error = 0;
x_n = x(1);
y_n = y(1);
if y(1) < y(2), ystep = 1; else ystep = -1; end 
myline = zeros(1,delx+1);
X = zeros(1,delx+1);
Y = zeros(1,delx+1);
for n = 1:delx+1
    if x_n <= 0 || y_n <= 0
        if n > 1
            myline(n) = myline(n-1);
            X(n) = X(n-1);
            Y(n) = Y(n-1);
        else
            myline(n) = 0;
            X(n) = 0;
            Y(n) = 0;
        end
    elseif steep
        myline(n) = mymat(y_n,x_n);
        %outmat(y_n,x_n) = outmat(y_n,x_n)+0.5;
        X(n) = y_n;
        Y(n) = x_n;
    else
        myline(n) = mymat(x_n,y_n);
        %outmat(x_n,y_n) = outmat(x_n,y_n)+0.5;
        X(n) = x_n;
        Y(n) = y_n;
    end    
    x_n = x_n + 1;
    error = error + dely;
    if 2*error >= delx, % same as -> if 2*error >= delx, 
        y_n = y_n + ystep;
        error = error - delx;
    end    
end
% -> a(y,x)
if dispFlag, imagesc(outmat); end
%plot(1:delx,myline)

function [q,r] = swap(s,t)
% function SWAP
q = t; r = s;
