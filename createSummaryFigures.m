function [ figHS ] = createSummaryFigures( configS )
%INITSUMMARYFIGURES Summary of this function goes here
%   Detailed explanation goes here

figHS = struct(...
    'summary',        [],...
    'map',            [],...
    'houghMask',      [],...
    'ellipseOutlier', [],...
    'planeOutlier',   []);

figHS(2).summary = [];

scrsz = get(0,'ScreenSize');
set(0,'Units','pixels');
% Outerposition: lower left corner X,Y, width, height

for i=1:2
    if configS.showSummaryFigsN >= 1 || configS.saveSummaryFigs == 1
        figHS(i).map     = figure('OuterPosition',...
            [0 scrsz(4)-scrsz(4)/1.5 scrsz(3)/2 scrsz(4)/1.5], 'Visible', 'off');
        figHS(i).summary = figure('OuterPosition',...
            [scrsz(3)/2 0 scrsz(3)/2 scrsz(4)], 'Visible', 'off');
    end
    
    if configS.showHoughMask == 1
        figHS(i).houghMask = figure('OuterPosition',...
            [scrsz(3)/4 0 scrsz(3)/4 scrsz(4)/3], 'Visible', 'off');
    end
    
    if configS.showEllipseOutlier == 1
        figHS(i).ellipseOutlier = figure('OuterPosition',...
            [0 0 scrsz(3)/4 scrsz(4)/3], 'Visible', 'off');
    end
    
    if configS.showPlaneOutlier == 1
        figHS(i).planeOutlier = figure('OuterPosition',...
            [0 0 scrsz(3)/3 scrsz(4)/3], 'Visible', 'off');
    end
 
end

end

