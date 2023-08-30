function [ configS ] = readConfigFile(dirName)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    %'pitFloorMeasuredDiameter', ...
names = {...
    'pitDiameter', ...
    'pitFloorDiameter', ...
    'pitFloorOverlapN', ...
    'polishedSurfaceInnerDiameter', ...
    'polishedSurfaceOuterDiameter',...
    'houghMaskInner',...
    'houghMaskOuter',...
    'ellipseOutlierMaxDist',...
    'ellipseOutlierDiscardNPoints',...
    'planeOutlierMaxDist',...
    'planeOutlierDiscardNPoints',...
    'showSummaryFigsN',...
    'showImgAnalFigs',...
    'showEllipseFigs',...
    'showPolishedElevFigs',...
    'showPitFloorFigs',...
    'isPolished',...
    'saveSummaryFigs',...
    'pitsNearGrainEdges',...
    'epoxyGrainElevationDifference',...
    'showPlaneOutlier', ...
    'firstEllipseFitting',...
    'secondEllipseFitting',...
    'showEllipseOutlier', ...
    'openUnprocessedZMP', ...
    'findRim', ...
    'keepPoints', ...
    'thetaB', ...
    'showHoughMask', ...
    'showPolishedSurfaceHistogram', ...
    'numTraining',...
    'showImgConvFigs',...
    'templateMaskOuter',...
    'templateMaskInner',...
    'findRimSigma',...
    'omitOutlierEllipsesSD',...
    'alpha',...
    'minimumFloorBlobSize',...
    'minimumBlobRelativeArea'%,...
    %'keepMBlobs'
    };

listing = dir([dirName filesep 'pitVolumeConfig*.txt']);
if length(listing) > 1
    error(['There is more than one configuration file in ' dirName]);
elseif isempty(listing)
    error(['Could not find a configuration file in ' dirName]);
end

fid = fopen([dirName filesep listing.name], 'r');
%fid = fopen([dirName filesep listing.name], 'r', 'n', 'UTF-8');
tline = fgetl(fid);
while ischar(tline)
    tline = strtrim(tline);
    % Notepad saves files with a UTF-BOM by default. That is a bit
    % annoying.
    %if(isempty(tline) || tline(1) == '%' || strncmp(tline, '﻿', 3))
    if(isempty(tline) || tline(1) == '%' || strncmp(tline, '﻿', 3))
    else
        C = textscan(tline, '%s = %f');
        
        
        if ~sum(strcmp(C{1}{1},names))
            error(['The configuration parameter "' C{1}{1} '" is not a '...
                'known parameter. Please change ' listing.name]);
        end
        
        configS.(C{1}{1}) = C{2};
    end
    
    tline = fgetl(fid);
end

fclose(fid);

%names = fieldnames(S);

ef = [dirName filesep listing.name];

for i=1:length(names)
    if ~isfield(configS, names{i})
        configError(ef, ['The parameter "' names{i} '" was not set.']);
    end
end

%% Basic error checking
%if configS.pitDiameter <= 0
%    configError(ef, '"pitDiameter" is less than 0.');
%end

if configS.pitFloorDiameter <= 0
    configError(ef, '"pitFloorDiameter" is less than 0.');
end

if configS.polishedSurfaceInnerDiameter <= 0
    configError(ef, '"polishedSurfaceInnerDiameter" is less than 0.');
end

if configS.polishedSurfaceOuterDiameter <= 0
    configError(ef, '"polishedSurfaceOuterDiameter" is less than 0.');
end

if configS.polishedSurfaceInnerDiameter > configS.polishedSurfaceOuterDiameter
    configError(ef, '"polishedSurfaceInnerDiameter" is larger than "polishedSurfaceOuterDiameter."');
end

if configS.houghMaskInner <= 0
    configError(ef, '"houghMaskInner" is less than 0.');
end

if configS.houghMaskOuter <= 0
    configError(ef, '"houghMaskOuter" is less than 0.');
end

if configS.houghMaskInner > configS.houghMaskOuter
    configError(ef, '"houghMaskInner" is larger than "houghMaskInner."');
end

if configS.ellipseOutlierMaxDist <= 0
    configError(ef, '"ellipseOutlierMaxDist" is less than 0.');
end

if configS.ellipseOutlierDiscardNPoints < 1
    configError(ef, '"ellipseOutlierDiscardNPoints" is less than 1. It must be an integer greater than 0.');
end

if configS.planeOutlierMaxDist <= 0
    configError(ef, '"planeOutlierMaxDist" is less than 0.');
end

if configS.planeOutlierMaxDist < 5
    disp(['"planeOutlierMaxDist" is quite small: ' num2str(configS.planeOutlierMaxDist) ' nm.']);
    disp('It should probably be >=5.');
    pause(3);
end

if configS.planeOutlierDiscardNPoints < 1
    configError(ef, '"planeOutlierDiscardNPoints" is less than 1. It must be an integer greater than 0.');
end

if configS.numTraining < 0
    configError(ef, '"numTraining" is less than 1. It should be 0 or >=2');
end
if configS.numTraining > 0 && configS.numTraining < 2
    configError(ef, ['"numTraining" is ' num2str(configS.numTraining) '. It should be 0 or >=2']);
end
if configS.pitFloorOverlapN > configS.numTraining
    configError(ef, ['"numTraining" is ' num2str(configS.numTraining) ...
        ' and "pitFloorOverlapN" is ' num2str(configS.pitFloorOverlapN) ...
        '. pitFloorOverlapN should be less than or equal to numTraining.']);
end


end

function configError(ef, str)
disp(['There is a problem with your configuration file: ' ef]);
error(str);
end
