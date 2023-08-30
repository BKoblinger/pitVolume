function [ pitGeomS, meanS, stdevS, rsdpS, meanOutS, stdevOutS, rsdpOutS, outlierIndex ] ...
    = getPitPopStatistics( pitGeomS, configS, dirName )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%{
pitGeomS = struct(...
    'label', '',...
    'filename', '',...
    'pitVol', '', ...
    'pitVolSD', '', ...
    'pitDepth', '', ...
    'pitDepthSD', '', ...
    'ellipseArea', '', ...
    'ellipseAreaSD', '', ...
    'polishedElev', '', ...
    'polishedElevSD', '', ...
    'polishedElevSDPlus', '', ...
    'polishedElevSDMinus', '', ...
    'pitFloorElev', '', ... 
    'pitFloorElevSD', '', ... 
    'ellipseA', '', ...
    'ellipseB', '', ...
    'ellipsePhi', '', ...
    'ellipseXc', '', ...
    'ellipseYc', '', ...
    'independentEllipse', '');
%}

meanS = nullPitGeomS();
stdevS = nullPitGeomS();
rsdpS = nullPitGeomS();
meanOutS = nullPitGeomS();
stdevOutS = nullPitGeomS();
rsdpOutS = nullPitGeomS();

X = [pitGeomS.('pitVol')];
nanIndex = isnan(X).*(1:length(X));
nanIndex = unique(nanIndex);
nanIndex(nanIndex == 0) = [];
if configS.saveSummaryFigs == 1 && ~isempty(nanIndex)
    copyDest = [dirName filesep 'ellipseFigs' filesep 'NaNs' filesep];
    mkdir(copyDest);
    disp('Copying NaNs ...');
    for i=1:length(nanIndex)
        [~,fileName,~] = fileparts(pitGeomS(nanIndex(i)).file);
        imageName = [dirName filesep 'ellipseFigs' filesep pitGeomS(nanIndex(i)).label '_' fileName];
        copySrc = [imageName '.tif'];
        fprintf(1, 'Source: %s; \tDestination: %s', copySrc, copyDest);
        [status,message] = copyfile(copySrc, copyDest);
        if status == 0
            disp('  ** Error! **  ');
        end
        fprintf(1,[message '\n']);
    end
end


    
SOutNames = {...
    'pitVol',...
    'pitDepth', ...
    'ellipseArea'};
outlierIndex = [];
for i = 1:numel(SOutNames)
    X = [pitGeomS.(SOutNames{i})];
    [outInd] = find_outliers_Thompson(X, configS.alpha);%, [], 1);
    outlierIndex = [outlierIndex outInd'];
    for j=1:length(outInd)
        pitGeomS(outInd(j)).outlier = [pitGeomS(outInd(j)).outlier SOutNames{i} ' '];
    end
end
outlierIndex = unique(outlierIndex);
if configS.saveSummaryFigs == 1 && ~isempty(outlierIndex)
    copyDest = [dirName filesep 'ellipseFigs' filesep 'outliers' filesep];
    mkdir(copyDest);
    disp('Copying outliers ...');
    for i=1:length(outlierIndex)
        [~,fileName,~] = fileparts(pitGeomS(outlierIndex(i)).file);
        imageName = [dirName filesep 'ellipseFigs' filesep pitGeomS(outlierIndex(i)).label '_' fileName];
        copySrc = [imageName '.tif'];
        fprintf(1, 'Source: %s; \tDestination: %s', copySrc, copyDest);
        [status,message] = copyfile(copySrc, copyDest);
        if status == 0
            disp('  ** Error! **  ');
        end
        fprintf(1,[message '\n']);
    end
end

X = [pitGeomS.('ellipseArea')];
X(isnan(X)) = [];
if configS.omitOutlierEllipsesSD == 1
    X(outlierIndex) = [];
end
ellipseAreaStd = std(X);
for i=1:length(pitGeomS)
    pitGeomS(i).ellipseAreaSD = ellipseAreaStd;
    pitGeomS(i).pitVolSD = pitGeomS(i).pitVol * sqrt((pitGeomS(i).pitDepthSD/pitGeomS(i).pitDepth)^2 ...
        + (pitGeomS(i).ellipseAreaSD/pitGeomS(i).ellipseArea)^2);
end

SNames = {...
    'pitVol',...
    'pitVolSD',...
    'pitDepth', ...
    'pitDepthSD', ...
    'ellipseArea', ...
    'ellipseAreaSD', ...
    'polishedElevSD', ...
    'polishedElevSDPlus', ...
    'polishedElevSDMinus', ...
    'pitFloorElevSD', ... 
    'ellipseA', ...
    'ellipseB', ...
    'ellipsePhi'};
for i = 1:numel(SNames)
    X = [pitGeomS.(SNames{i})];
    X(isnan(X)) = [];
    meanS.(SNames{i}) = mean(X);
    stdevS.(SNames{i}) = std(X);
    rsdpS.(SNames{i}) = 100*stdevS.(SNames{i}) / meanS.(SNames{i});

    X = [pitGeomS.(SNames{i})];
    X(outlierIndex) = [];
    X(isnan(X)) = [];
    meanOutS.(SNames{i}) = mean(X);
    stdevOutS.(SNames{i}) = std(X);
    rsdpOutS.(SNames{i}) = 100*stdevOutS.(SNames{i}) / meanOutS.(SNames{i});
end

meanS.label = 'Mean (with outliers)';
stdevS.label = '1 sd (with outliers)';
rsdpS.label = 'Re 1 sd (%) (with outliers)';
meanOutS.label = 'Mean';
stdevOutS.label = '1 sd';
rsdpOutS.label = 'Re 1 sd (%)';


end

