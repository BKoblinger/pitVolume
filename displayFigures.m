function [ ] = displayFigures( isVisible, label, stepN, figHS1, figHS2 )
%DISPLAYNEWFIGURES Summary of this function goes here
%   Detailed explanation goes here

% Play around with what figures will be displayed and closed
%figHS = struct(...
%    'summary', '',...
%    'map', '',...
%    'houghMask', '',...
%    'ellipseOutlier', '',...
%    'planeOutlier', '');

%if configS.showSummaryFigsN >= ind
if isVisible
    onoffStr = 'on';
else
    onoffStr = 'off';
end

set(0, 'CurrentFigure', figHS1.summary);
titleName = label;

% Add a backslash before the underscore in the string otherwise MATLAB
% treats it like a TeX command
k = strfind(titleName, '_');
if ~isempty(k)
    newTitle = [titleName(1:k(1)-1) '\'];
    for j=2:length(k)-1
        newTitle = [newTitle titleName(k(j):k(j+1)-1) '\'];
    end
    newTitle = [newTitle titleName(k(end):end)];
else
    newTitle = titleName;
end

suptitle(newTitle);

fields = fieldnames(figHS1);
for i=1:numel(fields)
    if ~((strcmp(fields{i}, 'planeOutlier') && (stepN == 0 || stepN == 1)) ...
        || (strcmp(fields{i}, 'houghMask') && (stepN == 1 || stepN == 2)) ...
        || (strcmp(fields{i}, 'ellipseOutlier') && (stepN == 2)))
    set(figHS1.(fields{i}), 'Name',label,'NumberTitle','off');    
    set(figHS1.(fields{i}), 'Visible', onoffStr);
    end
end

if nargin == 6 && ~isempty(figHS2)
    fields = fieldnames(figHS2);
    for i=1:numel(fields)
        set(figHS2.(fields{i}), 'Visible', 'off');
    end
end



%figHandles = findobj('Type','figure');
%if ~isempty(figHandles)
%    figHtoClose = intersect(figHandles, prevFigHandles);
%    prevFigHandles = setdiff(figHandles, prevFigHandles);
%    for j=1:length(figHtoClose)
%        close(figHtoClose(j));
%    end
%    %prevFigHandles = figHandles;
%end

end

