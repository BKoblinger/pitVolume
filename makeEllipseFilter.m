function [ ellipseFilter ] = makeEllipseFilter( mapS, pitGeomS, maskInner )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

xc = pitGeomS.ellipseXc;
yc = pitGeomS.ellipseYc;
phi = pitGeomS.ellipsePhi;
a = pitGeomS.ellipseA;
b = pitGeomS.ellipseB;

if nargin == 2
    maskInner = 0;
end

ellipseFilter = ...
    (((mapS.x-(xc))*cos(phi)+(mapS.y-(yc))*sin(phi)).^2/((a-maskInner)/mapS.scale)^2) + ...
    (((mapS.x-(xc))*sin(phi)-(mapS.y-(yc))*cos(phi)).^2/((b-maskInner)/mapS.scale)^2) <= 1;

end

