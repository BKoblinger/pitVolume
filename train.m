function [ corrMap, trainingMaps ] = train( dirName, listing, configS, figHSA )
%TRAIN Summary of this function goes here
%   Detailed explanation goes here


%a = zeros(1,configS.numTraining);
%b = zeros(1,configS.numTraining);
%phi = zeros(1,configS.numTraining);
cnt = 1;

trainingMaps{configS.numTraining} = [];
ellMap{configS.numTraining} = [];
indexes = 1:length(listing);
prevFigHandles = [];

figHS = figHSA(1);
%figHSprev = figHSA(2);

while cnt <= configS.numTraining
    if isempty(indexes)
        error('There are not enough images for the training.');
    end
    ind = indexes(1+round((length(indexes)-1)*rand()));
    
    %pitGeomS = nullPitGeomS();
    %pitGeomS.file = listing(ind).name;

    mapS = loadZMPfile([dirName filesep listing(ind).name], configS.openUnprocessedZMP);
    %pitGeomS.label = mapS.label;
    
    [ pitGeomS, ~, ~, ~, ~  ] ...
        = fitOneEllipse(mapS, [], mapS.label, listing(ind).name, configS, 1, figHS, []);
    
    %disp(['    Beginning ellipse training on: "' mapS.label '": ' listing(ind).name]);
    
    %plotOverviewMap( mapS.label, figHS, mapS, 1 );
    %
    %[ pitGeomS, boundary, borders, boundaryWithoutOutliers, isFit  ] = ...
    %    firstEllipseFitting( mapS, [], nullPitGeomS(), configS, figHS );

    %plotEllipseFit(borders, boundary, boundaryWithoutOutliers, pitGeomS, mapS, figHS);
    %plotPitEllipse(mapS, boundary, pitGeomS, mapS.scale, figHS);
    %displayFigures( 1, mapS.label, configS, 0, figHS );
    
    button = questdlg(['Is this a typical pit and is the ellipse well fit? ('...
        num2str(cnt) '/' num2str(configS.numTraining) ')']);
    if strcmp(button, 'Yes')
        trainingMaps{cnt} = mapS;
        ellMap{cnt}  = getEllipseMask(mapS, pitGeomS, configS);
        cnt = cnt+1;
    elseif strcmp(button, 'Cancel');
        error('You cancelled!');
    end
    indexes(indexes == ind) = [];
end

corrMap = convolveTrainingPits(ellMap, configS, mapS.scale);
%pause(5);

%close all;


end

