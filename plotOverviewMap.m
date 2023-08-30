function [  ] = plotOverviewMap( label, figHS, mapS, showOverviewMap )
%INITIALIZESUMMARYFIGURES Summary of this function goes here
%   Detailed explanation goes here

set(0,'CurrentFigure',figHS.map);
set(figHS.map, 'Name',label,'NumberTitle','off');
pitImageSc(mapS, 0, 'hot');

if showOverviewMap == 1
    set(figHS.map,'Visible','on');
end


end

