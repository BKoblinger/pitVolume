function pitVolume( varargin )
%%%%%
% A script to calculate the volume of pits formed by laser ablation and
% measured with the "ZeScope"
%
% General workflow:
%  - User input: directory with "digital elevation map" files and the 
%    configuration file
%  - A big loop (should be parallel processed; e.g. parfor)
%  - Import file (.zmp)
%  ** Find the area of the put
%  - Gaussian blur (remove high frequency noise)
%  - Edge detection
%  - Segmentation to get the edge of the pit as a matrix of (x,y) points
%  - Fit an ellipse to the data (least-squares method)
%  - Outlier detection
%  - Fit an ellipse to the outlier-removed edge
%  ** Get the average elevation of the polished surface
%  - ?? What is the appropriate radius (what are the characteristics of the
%  reprecipitated "junk"? Does it vary depending on pit size, depth, or
%  mineral? Probably best to get user input.
%  - Get the elevation in a series of "arcs"
%  - Elevation vs phi to a sine function (least squares)
%  - Outlier detection
%  - Refit sine function to less noisy data
%  ** Get the average elevation in the centre of the pit
%  - User input (in config file) for diameter
%  - Get average elevation in circle
%  ** Output
%  - Format requirements???
%%%%%%%

%% Initialize
%  - User input: directory with "digital elevation map" files and the 
%    configuration file
%  - A big loop (should be parallel processed; e.g. parfor)
%  - Import file (.zmp)



% Get the directory with all of the zmp files.
if isempty(varargin)
    % Ask for directory
    dirName = uigetdir();
    if dirName == 0
        error('You did not specify a directory.');
    end
else
    dirName = varargin{1};
end
%dirName = [dirName filesep 'zmp'];

% Read configuration file
 configS = readConfigFile(dirName);

% Get the names of all of the zmp files
listing = dir([dirName filesep '*.zmp']);
numFiles = length(listing);

if numFiles == 0
    error(['Error: No *.zmp files were found in the directory: "' dirName filesep '".  Exiting.']);
end

% Initialize the structure that will contains all of volume and geometry
% data
%pitGeomS(1) = nullPitGeomS();


figSummaryH = [];
subplotH = [];


if configS.saveSummaryFigs == 1
    if ~exist([dirName filesep 'ellipseFigs'],'dir')
        mkdir(dirName, 'ellipseFigs');
    end
    if ~exist([dirName filesep 'ellipseImages'], 'dir')
        mkdir(dirName, 'ellipseImages');
    end
end


%% First we find all of the ellipses then look for outliers!
[ pitGeomS, borders, boundary, boundaryWithoutOutliers, pitFloorModelS ]  = fitEllipsePopulation(dirName, listing, configS);
    
figHSA = createSummaryFigures(configS);
disp('STEP 3: FINDING ELEVATIONS');
for i=1:numFiles
    fprintf(1, ['Beginning elevation finding on file: "' pitGeomS(i).label '"; ' pitGeomS(i).file]);
    mapS = loadZMPfile([dirName filesep pitGeomS(i).file], configS.openUnprocessedZMP);
    %  ** Get the average elevation of the polished surface
    if mod(i,2) == 1
        figHS = figHSA(1);
        figHSprev = figHSA(2);
    else
        figHS = figHSA(2);
        figHSprev = figHSA(1);
    end
        
    if i <= configS.showSummaryFigsN
        plotOverviewMap( pitGeomS(i).label, figHS, mapS, 0 );
        if configS.saveSummaryFigs == 1
            set(0,'CurrentFigure',figHS.map);
            [~,fileName,~] = fileparts(pitGeomS(i).file);
            print(figHS.map, '-dtiff', [dirName filesep 'ellipseImages' filesep pitGeomS(i).label '_' fileName]);
            drawnow;
        end
    end
        
    if i <= configS.showSummaryFigsN || configS.saveSummaryFigs == 1
        plotEllipseFit(borders{i}, boundary{i}, boundaryWithoutOutliers{i}, ...
            pitGeomS(i), mapS, figHS);

        plotPitEllipse(mapS, boundary{i}, pitGeomS(i), mapS.scale, figHS);
    end
    
    fprintf(1, ' ... ');
      
    startTime = tic;
    try
        [pitGeomS(i).polishedElev, pitGeomS(i).polishedElevSD, pitGeomS(i).polishedElevSDPlus, pitGeomS(i).polishedElevSDMinus] = elevationPolishedSurface( ...
            pitGeomS(i).ellipseXc, pitGeomS(i).ellipseYc, configS, mapS, figHS);
    catch ME
        warning('Caught an error in elevationPolishedSurface');
        disp(ME);
        pitGeomS(i).polishedElev = nan;
        pitGeomS(i).polishedElevSD = nan;
        pitGeomS(i).polishedElevSDPlus = nan;
        pitGeomS(i).polishedElevSDMinus = nan;
    end
    pitGeomS(i).polishedElev = pitGeomS(i).polishedElev / 1000; % um
    pitGeomS(i).polishedElevSD = pitGeomS(i).polishedElevSD / 1000; % um
    pitGeomS(i).polishedElevSDPlus = pitGeomS(i).polishedElevSDPlus / 1000; % um
    pitGeomS(i).polishedElevSDMinus = pitGeomS(i).polishedElevSDMinus / 1000; % um
    ellapsedTime = toc(startTime);
    fprintf(1, '%.2f s\n', ellapsedTime);

    %  ** Get the average elevation in the centre of the pit
    % in nm
    %[pitGeomS(i).pitFloorElevOld, pitGeomS(i).pitFloorElevSDOld] = pitFloorElevationSimple( ...
    %    pitGeomS(i).ellipseXc, pitGeomS(i).ellipseYc, ...
    %    mapS, configS, figHS );
    try
        [pitGeomS(i).pitFloorElev, pitGeomS(i).pitFloorElevSD] = pitFloorElevation( ...
            mapS, configS, pitGeomS(i), pitFloorModelS, figHS );
    catch ME
        warning('Caught an error in pitFloorElevation');
        disp(ME);
        pitGeomS(i).pitFloorElev = nan;
        pitGeomS(i).pitFloorElevSD = nan;
    end
    %pitGeomS(i).pitFloorElevOld
    %pitGeomS(i).pitFloorElevSDOld
    %pitGeomS(i).pitFloorElev
    %pitGeomS(i).pitFloorElevSD

    pitGeomS(i).pitFloorElev = pitGeomS(i).pitFloorElev / 1000; % um
    pitGeomS(i).pitFloorElevSD = pitGeomS(i).pitFloorElevSD / 1000; % um

    pitGeomS(i).pitDepth = pitGeomS(i).polishedElev - pitGeomS(i).pitFloorElev; % um
    pitGeomS(i).pitDepthSD = sqrt(pitGeomS(i).pitFloorElevSD^2 + ...
        pitGeomS(i).polishedElevSD^2);

    % a / tan(thetaA) = b / tan(thetaB).  
    % V = (pi/3) a b (H^3 - (H-h)^3) / H^2, where H = b / tan(thetaB)

    if configS.thetaB == 0
        pitGeomS(i).pitVol = pitGeomS(i).pitDepth * pitGeomS(i).ellipseArea;
    else
        H = pitGeomS(i).ellipseB / tan(deg2rad(configS.thetaB));
        
        pitGeomS(i).pitVol = (pi/3) * pitGeomS(i).ellipseA * pitGeomS(i).ellipseB ...
            * (H^3 - (H - pitGeomS(i).pitDepth)^3) / H^2;
        %pitGeomS(i).pitVolSD = (pitGeomS(i).pitDepthSD/pitGeomS(i).pitDepth) * ...
            %pitGeomS(i).pitVol;
    end
    pitGeomS(i).pitVolSD = 0;
   

    if i <= configS.showSummaryFigsN || configS.saveSummaryFigs == 1
        displayFigures( i <= configS.showSummaryFigsN, pitGeomS(i).label, 2, figHS, figHSprev );
        %function [ ] = displayFigures( isVisible, label, configS, stepN, figHS1, figHS2 )
    end

    if configS.saveSummaryFigs == 1
        set(0, 'CurrentFigure', figHS.summary);
        [~,fileName,~] = fileparts(pitGeomS(i).file);
        print(figHS.summary, '-dtiff', [dirName filesep 'ellipseFigs' filesep pitGeomS(i).label '_' fileName]);
        drawnow;
    end

    %prevFigMapH = figMapH;
    %prevFigSummaryH = figSummaryH;
    if i == configS.showSummaryFigsN
        configS.showPlaneOutlier = 0;
        configS.showSummaryFigsN = 0;
    end
end

%figure();
%plot([pitGeomS.pitFloorElevSDOld]/1000, [pitGeomS.pitFloorElevSD], 'kx');
%axis square;
%axis equal;

[pitGeomS, meanS, stdevS, rsdpS, meanOutS, stdevOutS, rsdpOutS] ...
    = getPitPopStatistics(pitGeomS, configS, dirName);

%statsS = statsOutliers(pitGeomS);
%writePitVolumeFile(dirName, pitGeomS, statsS);
disp('Saving to file...');
writePitVolumeFile(dirName, pitGeomS, meanS, stdevS, rsdpS, meanOutS, stdevOutS, rsdpOutS);
disp('All finished!');
close all;



end
