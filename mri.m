%cd('/users/mgens/matlab_code/sharp-git');
info=dicominfo('mri_image.dcm');
voxelSize=zeros(1,2);
voxelSize(1:2)=info.PixelSpacing;
myImageUncropped=double(dicomread('mri_image.dcm'));
myMaskUncropped=double(dicomread('mri_mask.dcm'));

% display tumor
figure; imshow(myImageUncropped.*(myMaskUncropped+200.0),[]);
title('Slice from breast tumor MRI. Tumor in lighter shade.');

[row,col]=find(myMaskUncropped>0);
myImage=double(myImageUncropped(min(row):max(row),min(col):max(col)));
myMask=double(myMaskUncropped(min(row):max(row),min(col):max(col)));
[dimx dimy]=size(myImage);
sigmai_constant=2.5;
sigmax_constant=0.5;
whereTumor=find(myMask>0);
tumorIntensity=myImage(whereTumor);
sortedIntensity=sort(tumorIntensity);
quart1Intensity=sortedIntensity(ceil(size(sortedIntensity,1)*0.25));
quart3Intensity=sortedIntensity(ceil(size(sortedIntensity,1)*0.75));
sigmai=(quart3Intensity-quart1Intensity)*sigmai_constant/1.4; %interquartile / 1.4 = std dev for normal dist

% recursively partition tumor using Normalized Cut
[partsCell,partsCellLen,parentCell]=partition2D(myImage,myMask,sigmai,sigmax_constant,voxelSize);

% display first 8 levels of partitions
[whereTumorX,whereTumorY]=find(myMask>0);
myImageScaled=(myImage-min(tumorIntensity))/(max(tumorIntensity)-min(tumorIntensity));
myImageScaled=max(myImageScaled,0);
myImageScaled=min(myImageScaled,1);
numLevels=size(partsCell,1)-1;
figure('Position',[100 100 700 500]);
labelsTest2=zeros(size(whereTumorX,1),1);
currentLabel=1;
numRows=2;
myWidth=2;
for level=1:min(numLevels,8)
    for i=1:size(partsCell{level+1},1)/2
        voxels1=partsCell{level+1}((i*2)-1,1:partsCellLen{level+1}((i*2)-1));
        voxels2=partsCell{level+1}(i*2,1:partsCellLen{level+1}(i*2));
        labelsTest2(voxels1)=currentLabel;
        labelsTest2(voxels2)=currentLabel+1;
        currentLabel=currentLabel+2;
    end
    myCut=zeros(dimx,dimy);
    myCut(whereTumor)=labelsTest2;
    subplot(numRows,ceil(min(numLevels,8)/numRows),level);
    imshow(myImageScaled,[]);
    for i=1:dimx-1
        for j=1:dimy-1
            if myCut(i,j) ~= myCut(i+1,j)
                line([j-0.5 j+0.5],[i+0.5 i+0.5],'LineWidth',myWidth);
            end
            if myCut(i,j) ~= myCut(i,j+1)
                line([j+0.5 j+0.5],[i-0.5 i+0.5],'LineWidth',myWidth);
            end
        end
        if myCut(i,dimy) ~= myCut(i+1,dimy)
            line([dimy-0.5 dimy+0.5],[i+0.5 i+0.5],'LineWidth',myWidth);
        end
    end
    for j=1:dimy-1
        if myCut(dimx,j) ~= myCut(dimx,j+1)
            line([j+0.5 j+0.5],[dimx-0.5 dimx+0.5],'LineWidth',myWidth);
        end
    end
    title(gca,sprintf('%d',level));
end

% graph amount of heterogeneity versus scale
histBins=16;
sortedIntensity=sort(tumorIntensity);
minIntensity=sortedIntensity(ceil(size(sortedIntensity,1)*0.02));
maxIntensity=sortedIntensity(ceil(size(sortedIntensity,1)*0.98));
binCenters=minIntensity+[0.5:1:histBins-0.5]*(maxIntensity-minIntensity)/histBins;
emdList = [];

for level=1:numLevels
    for i=1:size(partsCell{level+1},1)/2
        voxels1=partsCell{level+1}((i*2)-1,1:partsCellLen{level+1}((i*2)-1));
        voxels2=partsCell{level+1}(i*2,1:partsCellLen{level+1}(i*2));
        counts1=hist(tumorIntensity(voxels1),binCenters);
        freqs1=counts1/sum(counts1);
        counts2=hist(tumorIntensity(voxels2),binCenters);
        freqs2=counts2/sum(counts2);
        emdList=[emdList; level,i,((size(voxels1,2)+size(voxels2,2))*prod(voxelSize))^(1/2), sum(abs(cumsum(freqs1)-cumsum(freqs2)))];
    end
end
numBinnedLevels=10;
x=zeros(numBinnedLevels,1); y=x;
for i=1:numBinnedLevels
    if i>1
        lowerBound=2^(i-3.5);
    else
        lowerBound=-10000;
    end
    if i<numBinnedLevels
        upperBound=2^(i-2.5);
    else
        upperBound=10000;
    end
    [ix,~]=find(squeeze(emdList(:,3)) >= lowerBound & squeeze(emdList(:,3)) < upperBound);
    if size(ix,1)
        x(i)=2^(i-3);
        y(i)=mean(emdList(ix,4));
    end
end
figure; semilogx(x,y,'b-x');
xlabel('Region scale (mm)');
ylabel('Signal intensity heterogeneity (mean EMD)');
title('Signal intensity heterogeneity as a function of scale');
