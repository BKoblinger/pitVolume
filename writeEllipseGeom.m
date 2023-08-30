function [  ] = writeEllipseGeom( dirName, pitGeomS, num )
%WRITEELLIPSEGEOM Summary of this function goes here
%   Detailed explanation goes here

fid = fopen([dirName filesep 'ellipseGeom_' num2str(num) '.tsv'], 'w');

sep = '\t';
fprintf(fid, [
    'Filename' sep ...
    'Ellipse area (um^2)' sep ...
    'Ellipse a (um)' sep ...
    'Ellipse b (um)' sep ...
    'Ellipse phi (rad)' sep ...
    'Ellipse x' sep ...
    'Ellipse y' sep ...
    'Ellipse fit independently\n']);

for i=1:length(pitGeomS)
    fprintf(fid, ['%s' sep], pitGeomS(i).file);
    fprintf(fid, ['%d' sep], pitGeomS(i).ellipseArea);
    fprintf(fid, ['%d' sep], pitGeomS(i).ellipseA);
    fprintf(fid, ['%d' sep], pitGeomS(i).ellipseB);
    fprintf(fid, ['%d' sep], pitGeomS(i).ellipsePhi);
    fprintf(fid, ['%d' sep], pitGeomS(i).ellipseXc);
    fprintf(fid, ['%d' sep], pitGeomS(i).ellipseYc);
    if ~isempty(pitGeomS(i).independentEllipse)
        fprintf(fid, '%s', pitGeomS(i).independentEllipse);
    end
    fprintf(fid, '\n');
end

fclose(fid); 


end

