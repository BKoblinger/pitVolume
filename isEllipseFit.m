function [ isFit ] = isEllipseFit( object, msg )
%COULDNOTFITELLIPSE Summary of this function goes here
%   Detailed explanation goes here
isFit = true;
if isempty(object) || (isstruct(object) && isnan(object.a))
    isFit = false;
    disp(['    *** Could not fit an ellipse. ' msg]);

end

end

