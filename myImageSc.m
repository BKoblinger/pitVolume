function [  ] = myImageSc( Map )
%MYIMAGESC Summary of this function goes here
%   Detailed explanation goes here

%use the following to clip the displayed Z-range from -1%PV to +1%PV
%colorlimiter = [-0.01*PV 0.01*PV];
%disable the following, empty colorlimiter, when using line above
%colorlimiter = [minMap maxMap];
%prepare map display in figure window
%therefore set units to pixels
%set(0,'Units','pixels');
%scrsz = get(0,'ScreenSize');
%generate a figure window in the middle of the screen of half size
%of the screen
%figure('Position',[scrsz(3)/4 scrsz(4)/4 scrsz(3)/2 scrsz(4)/2]);
figure();
%display map
%imagesc(Map,colorlimiter);
imagesc(Map);
set(gca, 'YDir', 'normal');
axis equal;
%title({['Map file: ' filename], ...
%    ['Min: ' num2str(minMap,'%7.3f') ' nm, Max: ' num2str(maxMap,'%7.3f') ' nm, PV: ' num2str(PV,'%7.3f') ' nm, RMS: ' num2str(RMS,'%7.3f') ' nm'], ...
%    ['Number of void pixels: ' num2str(voidPixels) ' , number of valid Pixels: ' num2str(validPixels)]});
colormap(hot);
hcb = colorbar();
set(get(hcb,'Title'),'String','Elevation (nm)');
drawnow();


end

