function [ pitGeomS, borders, boundary, boundaryWithoutOutliers, pitFloorModelS ] ...
    = fitEllipsePopulation( dirName, listing, configS )
%FITELLIPSEPOPULATION Summary of this function goes here
%   Detailed explanation goes here

%prevFigHandles = [];
boundary = cell(length(listing),1);
borders = cell(length(listing),1);
boundaryWithoutOutliers = cell(length(listing),1);

figHSA = createSummaryFigures( configS );

%% TRAINING
isNewTraining = questdlg('Would you like to have a new model pit training session, or load an existing one?','New training session?','New','Load','Cancel','New');
switch isNewTraining
    case 'New'
        %% Determine the average ellipse mask
        convMap = [];
        trainingMaps = [];
        if configS.numTraining >= 1 && length(listing) >= configS.numTraining ...
                && (configS.firstEllipseFitting == 1 || configS.secondEllipseFitting == 1)
            disp('STEP 1A: TRAINING');
            [convMap, trainingMaps] = train(dirName, listing, configS, figHSA);

        elseif configS.pitDiameter == 0
            error(['You specified ' num2str(configS.numTraining) ...
                ' ellipse were needed (numTraining), but there is/are only '...
                num2str(length(listing)) ' .zmp file(s) in directory: ' dirName]);
        end

        %% Make a model of the pit floor
        if ~isempty(trainingMaps)
            disp('STEP 1B: PIT FLOOR MODEL');
            pitFloorModelS = makeAvgPitFloor(trainingMaps, convMap, configS, figHSA);
            uisave({'convMap','pitFloorModelS'}, [dirName filesep 'modelPit.mat']);
        else
            disp('    Cannot make a model of the pit floor because there are no training ellipses.');
        end
    case 'Load'
        uiopen([dirName filesep 'modelPit.mat']);
    case 'Cancel'
        exit(1);
end

%% PROCESSING
%% Now fit all the ellipses!
pitGeomS(length(listing)) = nullPitGeomS();
for i = 1:length(listing)
    pitGeomS(i) = nullPitGeomS();
end

if configS.firstEllipseFitting == 1
    disp('STEP 2: ELLIPSE FITTING');
    for i=1:length(listing)
    %parfor i = 1:numFile
        if mod(i,2) == 1
            figHS = figHSA(1);
            figHSprev = figHSA(2);
        else
            figHS = figHSA(2);
            figHSprev = figHSA(1);
        end

        %% Load the file, etc
        mapS = loadZMPfile([dirName filesep listing(i).name],...
            configS.openUnprocessedZMP);

        try
            [ pitGeomS(i), boundary{i}, borders{i}, boundaryWithoutOutliers{i}, isFit  ] = ...
                fitOneEllipse( mapS, convMap, mapS.label, listing(i).name, configS, i <= configS.showSummaryFigsN, figHS, figHSprev );
        catch ME
            warning('Caught an error in fitOneEllipse');
            disp(ME);
            pitGeomS(i) = nan;
            boundary{i} = nan;
            borders{i} = nan;
            boundaryWithoutOutliers{i} = nan;
            isFit = 0;
        end
        
        if i == configS.showSummaryFigsN
            configS.showHoughMask = 0;
            configS.showEllipseOutlier = 0;
            close all;
            drawnow;
        end

    end
    
    save([dirName filesep 'ellipseGeom1.mat'],'pitGeomS','boundary','boundaryWithoutOutliers', 'borders', 'pitFloorModelS');
    %writeEllipseGeom(dirName, pitGeomS, 1);
else
   load([dirName filesep 'ellipseGeom1.mat'],'pitGeomS','boundary','boundaryWithoutOutliers', 'borders', 'pitFloorModelS');
   if isempty(boundary) || isempty(boundaryWithoutOutliers)
       error(['Did not find "boundary" or "boundaryWithoutOutliers" in ellipseGeom1.mat. ' ...
           'Try re-running but first set the parameter firstEllipseFitting = 1.']);
   end
end
disp('STEP 2 FINISHED.');

[outlierIndex, meanA, stdA, meanB, stdB, meanPhi] = findOutlierEllipses(pitGeomS);

outlierInd = find(outlierIndex);

if isnan(meanA) || isempty(meanA)
    error('***** Could not find ellipses in any of the files. Sorry. *****');
end

%i = [];
if configS.secondEllipseFitting == 1 
    % Loop over each file
    disp('STEP 2B: SECOND ELLIPSE FITTING AND ELEVATION FINDING');
    disp('**** WARNING: THE SECOND ELLIPSE FITTING IS CURRENTLY NOT RECOMMENDED. ****');
    disp('**** TALK TO BRETT IF YOU WANT IT. ****');
    disp('**** IT WILL PROBABLY CRASH ... ***');
    pause(60);
    
    configS.findRim = 0;
        %parfor i = 1:numFile
    for j = 1:length(outlierInd)
        ind = outlierInd(j);
        %[~,i] = max(outliers);
        [ figSummaryH, figMapH ] = initSummaryFigures( pitGeomS(i).label, configS, mapS, 1 );
        tic
        disp(['    Beginning second ellipse fitting calculations on file: ' pitGeomS(ind).filename]);

        mapS = loadZMPfile([dirName filesep pitGeomS(ind).filename], configS.openUnprocessedZMP);


        if j <= configS.showSummaryFigsN
            set(0,'CurrentFigure',figMapH);
            pitImageSc(mapS, 0, 'hot');
        end

        tic
        %  ** Find the area of the pit
        %  - Segmentation to get the edge of the pit as a matrix of (x,y) points
        pitGeomS(ind).filename

        [xc, yc, boundary{ind}] = ...
            fitASpecifiedEllipse(mapS, ...
            pitGeomS(ind).ellipseXc, pitGeomS(ind).ellipseYc, meanA, meanB, ...
            meanPhi, convMap, configS);
        %[pitGeomS(i).ellipseXc, pitGeomS(i).ellipseYc, boundary{i}] = ...
            %fitASpecifiedEllipse(mapS, boundary{i}, pitGeomS(i).ellipseXc, pitGeomS(i).ellipseYc, medianA, medianB, medianPhi);
        %boundary{i} = pitSegmentation(mapS, configS);
        %pitGeomS(i).ellipseA = medianA; % in um
        %pitGeomS(i).ellipseB = medianB; % in um
        %pitGeomS(i).ellipsePhi = medianPhi; % in rad
        %pitGeomS(i).ellipseArea = pi * medianA * medianB;  % in um^2


        %if configS.showSummaryFigsN == 1
        %    set(0,'CurrentFigure',figSummaryH);
        %    subplotH = subplot(3,2,2,'replace');
        %end

        % Find the ellipse area in pixels and the ellipse parameters
        [ellipseS, boundaryWithoutOutliers{ind}] = pitArea(boundary{ind}, configS, mapS.scale, 1);    
        pitGeomS(ind).ellipseA = ellipseS.a * (mapS.scale); % in um
        pitGeomS(ind).ellipseB = ellipseS.b * (mapS.scale); % in um
        pitGeomS(ind).ellipsePhi = ellipseS.phi; % in rad
        pitGeomS(ind).ellipseXc = ellipseS.X0_in; % in pixels!
        pitGeomS(ind).ellipseYc = ellipseS.Y0_in; % in pixels!
        pitGeomS(ind).ellipseArea = pi * ellipseS.a * ellipseS.b;  % in um^2
        toc

        if abs(pitGeomS(ind).ellipseA - meanA) > 2.5*stdA ...
            || abs(pitGeomS(ind).ellipseB - meanB) > 2.5*stdB
            disp('    Second ellipse fitting was unsuccessful. Setting to median values.');
            pitGeomS(ind).ellipseA = meanA; % in um
            pitGeomS(ind).ellipseB = meanB; % in um
            pitGeomS(ind).ellipsePhi = meanPhi; % in rad
            pitGeomS(ind).ellipseXc = xc; % in pixels!
            pitGeomS(ind).ellipseYc = yc; % in pixels!
            pitGeomS(ind).ellipseArea = pi * meanA * meanB;  % in um^2
            pitGeomS(ind).independentEllipse = 'FALSE';
        end

        if j <= configS.showSummaryFigsN
            set(0,'CurrentFigure',figSummaryH);
            subplotH = subplot(3,2,2,'replace');
            plotEllipseFit(borders{ind}, boundary{ind}, ...
                boundaryWithoutOutliers{ind}, pitGeomS(ind), mapS, ...
                configS);

            subplotH = subplot(3,2,1,'replace');
            plotPitEllipse(subplotH, mapS, boundary{ind}, pitGeomS(ind), ...
                mapS.scale);
        end
        
         displayFigures( configS, 1, figHS, figHSprev );
         toc

    end
    save([dirName filesep 'ellipseGeom2.mat'],'pitGeomS','boundary','boundaryWithoutOutliers', 'borders');
    writeEllipseGeom(dirName, pitGeomS, 2);
   
%else
%    load([dirName filesep 'ellipseGeom2.mat'],'pitGeomS','boundary','boundaryWithoutOutliers', 'borders');
%    if isempty(boundary) || isempty(boundaryWithoutOutliers)
%        error(['Did not find "boundary" or "boundaryWithoutOutliers" in ellipseGeom2.mat. ' ...
%           'Try re-running but first set the parameter firstEllipseFitting = 1.']);
%    end
end
disp('STEP 2 FINISHED');
close all;


