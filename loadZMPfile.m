function [ MapS ] = loadZMPfile( filename, openUnprocessedZMP )
%loadZMPfile Load and return Zemetric .zmp
%   Provides file transfer functions for Matlab to read Zemetrics
%   .zmp-files.
%   Returns:
%     Map - a matrix of elevations in nm
%     scale - Spatial resolution, mm/pix

% ReadZMP.m
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% provides file transfer functions for Matlab to read Zemetrics .zmp-files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2012/04/02 Wolfgang Kaehler
% (c) 2012 by Zygo/Zemetrics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%clear all;            % make workspace empty
%clc;                  % make screen empty
%close all;            % close all windows with figures  
%warning off MATLAB:divideByZero
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%declaration/defines



verbosity = false;



%{
[filename,filepath] = uigetfile('*.zmp','Please select ZMP file',filepath);
if ~isequal(filename, 0)
    filename=strcat(filepath,filename);
    disp(filename);
else
    return;
end
%}



%disable the previous lines with the file dialog and enable the following
%to load the test map file in the same directory as this source file

%filename = 'test.zmp';

MLVer=ver('Matlab');
[mainver,remain]=strtok(MLVer.Version,'.');
subver=strtok(remain,'.');
if (str2double(mainver) >= 7) && (str2double(subver) > 5)
    [lError, InfoHeader, Map]=ReadZMPFileME(filename, verbosity);
else
    [lError, InfoHeader, Map, sequenceLabel]=ReadZMPFile(filename, verbosity, openUnprocessedZMP);
end

if lError==0
    %disp('file read successfully.');
    %file successfully loaded, display information
    %from this point on you can write your own code to get information
    %from the map, display it or parts of it
    %disp('InfoHeader =');
    %disp(InfoHeader);
    %disp (['Spatial resolution in X = ' num2str(InfoHeader.XScale) ' mm/pix']);
    %disp (['Spatial resolution in Y = ' num2str(InfoHeader.YScale) ' mm/pix']);
    if InfoHeader.XScale == InfoHeader.YScale
        scale = InfoHeader.XScale;
    else
        disp(['X and Y scale are not the same: ' filename]);
    end
    
    %calculate map size x,y in mm
    %mapSizeX = double(InfoHeader.Cols) * InfoHeader.XScale;
    %mapSizeY = double(InfoHeader.Rows) * InfoHeader.YScale;
    %disp (['Field size in X = ' num2str(mapSizeX) ' mm']);
    %disp (['Field size in Y = ' num2str(mapSizeY) ' mm']);
    
    %calculate RMS from valid pixels, void pixels are NaN
    %RMS=std(Map(~isnan(Map(:))));
    %calculate minimum
    %minMap = min(Map(~isnan(Map(:))));
    %calculate maximum
    %maxMap = max(Map(~isnan(Map(:))));
    %calculate peak-to-valley
    %PV=maxMap-minMap;

    %count void pixels
    %voidPixels = sum(isnan(Map(:)));
    %count valid pixels
    %validPixels = sum(~isnan(Map(:)));
    %disp(['Number of void pixels: ' num2str(voidPixels) ' , number of valid Pixels: ' num2str(validPixels)]);
    
    %use the following to clip the displayed Z-range from -1%PV to +1%PV
    %%colorlimiter = [-0.01*PV 0.01*PV];
    %%disable the following, empty colorlimiter, when using line above
    %colorlimiter = [minMap maxMap];
    %%prepare map display in figure window
    %%therefore set units to pixels
    %set(0,'Units','pixels');
    %scrsz = get(0,'ScreenSize');
    %%generate a figure window in the middle of the screen of half size
    %%of the screen
    %figure('Position',[scrsz(3)/4 scrsz(4)/4 scrsz(3)/2 scrsz(4)/2]);
    %%display map
    %imagesc(Map,colorlimiter);
    %axis equal;
    %title({['Map file: ' filename], ...
    %    ['Min: ' num2str(minMap,'%7.3f') ' nm, Max: ' num2str(maxMap,'%7.3f') ' nm, PV: ' num2str(PV,'%7.3f') ' nm, RMS: ' num2str(RMS,'%7.3f') ' nm'], ...
    %    ['Number of void pixels: ' num2str(voidPixels) ' , number of valid Pixels: ' num2str(validPixels)]});
    %colormap(hot);
    %colorbar;
else
    disp(['error reading file : ' filename]);
end

if isempty(sequenceLabel)
    [~,fname,~] = fileparts(filename);
    sequenceLabel = fname;
end

MapS.label = sequenceLabel;
[colsInImage, rowsInImage] = meshgrid(1:length(Map(1,:)), 1:length(Map(:,1)));
MapS.mapZ = Map;
MapS.scale = scale * 1000;
MapS.x = colsInImage;
MapS.y = rowsInImage;


end

