function [ ] = writePitVolumeFile( dirName, pitGeomS, meanS, stdevS, rsdpS, meanOutS, stdevOutS, rsdpOutS )
%WRITEPITVOLUMEFILE Summary of this function goes here
%   Detailed explanation goes here

fid = fopen([dirName filesep 'pitVolumesOut.tsv'], 'w');

sep = '\t';
fprintf(fid, [
    'Label' sep ...
    'Filename' sep ...
    'Outlier' sep ...
    'Pit volume (um^3)' sep ...
    'Pit volume 1 sd (um^3)' sep ...
    'Pit volume rsd (%%)' sep ...
    'Pit depth (um)' sep ...
    'Pit depth 1 sd (um)' sep ...
    'Pit depth rsd (%%)' sep ...
    'Ellipse area (um^2)' sep ...
    'Ellipse area 1 sd (um^2)' sep ...
    'Ellipse area rsd (%%)' sep ...
    'Polished surface elevation (um)' sep ...
    'Polished surface elevation 1 sd (um)' sep ...
    'Polished surface elevation 1 sd positive (um)' sep ...
    'Polished surface elevation 1 sd negative (um)' sep ...
    'Pit floor elevation (um)' sep ...
    'Pit floor elevation 1 sd (um)' sep ...
    'Ellipse a (um)' sep ...
    'Ellipse b (um)' sep ...
    'Ellipse phi (rad)\n']);

pitGeomS(end+1) = meanOutS;
pitGeomS(end+1) = stdevOutS;
pitGeomS(end+1) = rsdpOutS;
pitGeomS(end+1) = meanS;
pitGeomS(end+1) = stdevS;
pitGeomS(end+1) = rsdpS;

for i=1:length(pitGeomS)
    if ~isempty(strfind(pitGeomS(i).label, 'Mean'))
        fprintf(fid, '\n');
    end
    fprintf(fid, ['%s' sep], pitGeomS(i).label);
    fprintf(fid, ['%s' sep], pitGeomS(i).file);
    fprintf(fid, ['%s' sep], pitGeomS(i).outlier);
    fprintf(fid, ['%d' sep], pitGeomS(i).pitVol);
    fprintf(fid, ['%d' sep], pitGeomS(i).pitVolSD);
    fprintf(fid, ['%.3f' sep], 100*pitGeomS(i).pitVolSD/pitGeomS(i).pitVol);
    fprintf(fid, ['%d' sep], pitGeomS(i).pitDepth);
    fprintf(fid, ['%d' sep], pitGeomS(i).pitDepthSD);
    fprintf(fid, ['%.3f' sep], 100*pitGeomS(i).pitDepthSD/pitGeomS(i).pitDepth);
    fprintf(fid, ['%d' sep], pitGeomS(i).ellipseArea);
    fprintf(fid, ['%d' sep], pitGeomS(i).ellipseAreaSD);
    fprintf(fid, ['%.3f' sep], 100*pitGeomS(i).ellipseAreaSD/pitGeomS(i).ellipseArea);
    fprintf(fid, ['%d' sep], pitGeomS(i).polishedElev);
    fprintf(fid, ['%d' sep], pitGeomS(i).polishedElevSD);
    fprintf(fid, ['%d' sep], pitGeomS(i).polishedElevSDPlus);
    fprintf(fid, ['%d' sep], pitGeomS(i).polishedElevSDMinus);
    fprintf(fid, ['%d' sep], pitGeomS(i).pitFloorElev);
    fprintf(fid, ['%d' sep], pitGeomS(i).pitFloorElevSD);
    fprintf(fid, ['%d' sep], pitGeomS(i).ellipseA);
    fprintf(fid, ['%d' sep], pitGeomS(i).ellipseB);
    fprintf(fid, ['%d' sep], pitGeomS(i).ellipsePhi);
    if ~isempty(pitGeomS(i).independentEllipse)
        fprintf(fid, '%s', pitGeomS(i).independentEllipse);
    end
    fprintf(fid, '\n');
end

fclose(fid); 


end

