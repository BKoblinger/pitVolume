function [ corrMap ] = convolveTrainingPits( maps, configS, scale )
%CORRELATEIMAGES Summary of this function goes here
%   Detailed explanation goes here

center = round(size(maps{1})/2);
xc = center(2);
yc = center(1);

trimmed{length(maps)} = [];
for i=1:length(maps)
    trimmed{i} = maps{i};
    J = trimmed{i} == 0;
    trimmed{i}(:,all(J,1))='';
    trimmed{i}(all(J,2),:)='';
end

x = zeros(length(maps));
y = zeros(length(maps));
for i=1:length(maps)
    for j=1:length(maps)
        if i==j
            continue;
        end

        C12 = conv2(maps{i},trimmed{j}, 'same');

        [~, ind] = max(C12(:));
        [y(i,j), x(i,j)] = ind2sub(size(C12),ind);

        %if configS.showImgConvFigs == 1
        %    myImageSc(maps{i});
        %    myImageSc(trimmed{j});
        %    myImageSc(C12);
        %end
    end
end

shifted{length(maps)} = [];
for i=1:length(maps)
    xAvg = x(i,:);
    xAvg(xAvg == 0) = [];
    xAvg = mean(xAvg);
    yAvg = y(i,:);
    yAvg(yAvg == 0) = [];
    yAvg = mean(yAvg);
    shifted{i} = mytranslate(maps{i},[xc-xAvg, yc-yAvg]);
    
    if configS.showImgConvFigs == 1
        myImageSc(maps{i});
        myImageSc(shifted{i});
    end
end
  
corrMap = shifted{1};
for i=2:length(shifted)
    %corrMap = corrMap & shifted{i};
    corrMap = corrMap + shifted{i};
end
corrMap = corrMap/length(shifted);
%myImageSc(corrMap);
fillCorrMap = imfill(1.0*corrMap);

SE = strel('disk', round(configS.houghMaskOuter/scale));
rimMask = imdilate(fillCorrMap, SE);
rimMask = rimMask - fillCorrMap;
mask = 3.0*corrMap - 1.0*fillCorrMap - 1.0*rimMask;
%myImageSc(mask)


J = mask == 0;
mask(:,all(J,1))='';
mask(all(J,2),:)='';
if configS.showImgConvFigs == 1
    myImageSc(mask);
    pause(5);
end

corrMap = mask;
%myImageSc(corrMap);

end

