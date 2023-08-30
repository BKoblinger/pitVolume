function [ ellMap ] = getEllipseMask( mapS, pitGeomS, configS )
%GETELLIPSEMASK Summary of this function goes here
%   Detailed explanation goes here

ellipseFilter = makeEllipseFilter(mapS, pitGeomS);

nanMap = isnan(mapS.mapZ);
ellMap = ellipseFilter .* nanMap;

if configS.showImgAnalFigs == 1
    myImageSc(ellipseFilter);
    myImageSc(nanMap);
    myImageSc(ellMap);
end

end

