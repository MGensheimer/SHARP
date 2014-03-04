info=dicominfo('mri_image.dcm');
voxelSize=zeros(1,2);
voxelSize(1:2)=info.PixelSpacing;
imageUncropped=double(dicomread('mri_image.dcm'));
maskUncropped=double(dicomread('mri_mask.dcm'));

% display tumor
figure; imshow(imageUncropped.*(maskUncropped+200.0),[]);
title('Slice from breast tumor MRI. Tumor in lighter shade.');

[row,col]=find(maskUncropped>0);
image=double(imageUncropped(min(row):max(row),min(col):max(col)));
mask=double(maskUncropped(min(row):max(row),min(col):max(col)));
mask(mask<0)=0;
[dimx dimy]=size(image);
sigmai_constant=2.5;
sigmax_constant=0.5;
whereTumor=find(mask>0);
[whereTumorX,whereTumorY]=find(mask>0);
posToIndex=zeros(dimx,dimy);
posToIndex(whereTumor)=1:size(whereTumor,1);
tumorIntensity=image(whereTumor);
sortedIntensity=sort(tumorIntensity);
quart1Intensity=sortedIntensity(ceil(size(sortedIntensity,1)*0.25));
quart3Intensity=sortedIntensity(ceil(size(sortedIntensity,1)*0.75));
sigmai=(quart3Intensity-quart1Intensity)*sigmai_constant/1.4; %interquartile / 1.4 = std dev for normal dist



% recursively partition tumor using Normalized Cut
[partsCell,partsCellLen,parentCell]=partition2D(image,mask,sigmai,sigmax_constant,voxelSize);



% display first 8 levels of partitions
padding=10;
imagePadded=double(imageUncropped(min(row)-padding:max(row)+padding,min(col)-padding:max(col)+padding));
imageScaled=(imagePadded-min(tumorIntensity))/(max(tumorIntensity)-min(tumorIntensity));
imageScaled=max(imageScaled,0);
imageScaled=min(imageScaled,1);
numLevels=size(partsCell,1)-1;
figure('Position',[100 100 700 500]);
labelsList=zeros(size(whereTumorX,1),1);
currentLabel=1;
numRows=2;
myWidth=2;
for level=1:min(numLevels,8)
    for i=1:size(partsCell{level+1},1)/2
        voxels1=partsCell{level+1}((i*2)-1,1:partsCellLen{level+1}((i*2)-1));
        voxels2=partsCell{level+1}(i*2,1:partsCellLen{level+1}(i*2));
        labelsList(voxels1)=currentLabel;
        labelsList(voxels2)=currentLabel+1;
        currentLabel=currentLabel+2;
    end
    cut=zeros(dimx,dimy);
    cut(whereTumor)=labelsList;
    subplot(numRows,ceil(min(numLevels,8)/numRows),level);
    imshow(imageScaled,[]);
    for i=1:dimx-1
        for j=1:dimy-1
            if cut(i,j) ~= cut(i+1,j)
                line([j-0.5+padding j+0.5+padding],[i+0.5+padding i+0.5+padding],'LineWidth',myWidth);
            end
            if cut(i,j) ~= cut(i,j+1)
                line([j+0.5+padding j+0.5+padding],[i-0.5+padding i+0.5+padding],'LineWidth',myWidth);
            end
        end
        if cut(i,dimy) ~= cut(i+1,dimy)
            line([dimy-0.5+padding dimy+0.5+padding],[i+0.5+padding i+0.5+padding],'LineWidth',myWidth);
        end
    end
    for j=1:dimy-1
        if cut(dimx,j) ~= cut(dimx,j+1)
            line([j+0.5+padding j+0.5+padding],[dimx-0.5+padding dimx+0.5+padding],'LineWidth',myWidth);
        end
    end
    for i=1:dimx %fill in edges
        if cut(i,1)
            line([0.5 0.5]+padding,[i-0.5 i+0.5]+padding,'LineWidth',myWidth);
        end
        if cut(i,dimy)
            line([dimy+0.5 dimy+0.5]+padding,[i-0.5 i+0.5]+padding,'LineWidth',myWidth);
        end                            
    end
    for j=1:dimy
        if cut(1,j)
            line([j-0.5 j+0.5]+padding,[0.5 0.5]+padding,'LineWidth',myWidth);
        end
        if cut(dimx,j)
            line([j-0.5 j+0.5]+padding,[dimx+0.5 dimx+0.5]+padding,'LineWidth',myWidth);
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



%find heterogeneity versus scale at a specified point
queryPoint=[18,28];
emdList = [];
labelsList=zeros(size(whereTumorX,1),1);
currentLabel=1;
numLevels=size(partsCell,1)-1;

padding=10;
imagePadded=double(imageUncropped(min(row)-padding:max(row)+padding,min(col)-padding:max(col)+padding));
imageScaled=(imagePadded-min(tumorIntensity))/(max(tumorIntensity)-min(tumorIntensity));
imageScaled=max(imageScaled,0);
imageScaled=min(imageScaled,1);
figure('Position',[100 100 700 500]);
%cut=mask;
imshow(imageScaled,[]);
myWidth=2;
colorList=lines(numLevels+1);
for level=1:numLevels
    for i=1:size(partsCell{level+1},1)/2
        voxels1=partsCell{level+1}((i*2)-1,1:partsCellLen{level+1}((i*2)-1));
        voxels2=partsCell{level+1}(i*2,1:partsCellLen{level+1}(i*2));
        if max(voxels1==posToIndex(queryPoint(1),queryPoint(2))) || max(voxels2==posToIndex(queryPoint(1),queryPoint(2)))
            counts1=hist(tumorIntensity(voxels1),binCenters);
            freqs1=counts1/sum(counts1);
            counts2=hist(tumorIntensity(voxels2),binCenters);
            freqs2=counts2/sum(counts2);
            emdList=[emdList; level,i,((size(voxels1,2)+size(voxels2,2))*prod(voxelSize))^(1/2), sum(abs(cumsum(freqs1)-cumsum(freqs2)))];

            labelsList(voxels1)=currentLabel;
            labelsList(voxels2)=currentLabel+1;
            currentLabel=currentLabel+2;

            cut=zeros(dimx,dimy); cutEroded=cut;
            cut(whereTumor(voxels1))=1;
            cut(whereTumor(voxels2))=2;
            for i=1:dimx-1
                for j=1:dimy-1
                    if cut(i,j) + cut(i+1,j)==3
                        line([j-0.5+padding j+0.5+padding],[i+0.5+padding i+0.5+padding],'LineWidth',myWidth,'Color',squeeze(colorList(level+1,:)));
                    end
                    if cut(i,j) + cut(i,j+1)==3
                        line([j+0.5+padding j+0.5+padding],[i-0.5+padding i+0.5+padding],'LineWidth',myWidth,'Color',squeeze(colorList(level+1,:)));
                    end
                end
                if cut(i,dimy) + cut(i+1,dimy)==3
                    line([dimy-0.5+padding dimy+0.5+padding],[i+0.5+padding i+0.5+padding],'LineWidth',myWidth,'Color',squeeze(colorList(level+1,:)));
                end
            end
            for j=1:dimy-1
                if cut(dimx,j) + cut(dimx,j+1)==3
                    line([j+0.5+padding j+0.5+padding],[dimx-0.5+padding dimx+0.5+padding],'LineWidth',myWidth,'Color',squeeze(colorList(level+1,:)));
                end
            end
        end
    end
end
% show outline of whole tumor
cut=mask;
for i=1:dimx-1
    for j=1:dimy-1
        if cut(i,j) ~= cut(i+1,j)
            line([j-0.5+padding j+0.5+padding],[i+0.5+padding i+0.5+padding],'LineWidth',myWidth,'Color',squeeze(colorList(1,:)));
        end
        if cut(i,j) ~= cut(i,j+1)
            line([j+0.5+padding j+0.5+padding],[i-0.5+padding i+0.5+padding],'LineWidth',myWidth,'Color',squeeze(colorList(level+1,:)));
        end
    end
    if cut(i,dimy) ~= cut(i+1,dimy)
        line([dimy-0.5+padding dimy+0.5+padding],[i+0.5+padding i+0.5+padding],'LineWidth',myWidth,'Color',squeeze(colorList(1,:)));
    end
end
for j=1:dimy-1
    if cut(dimx,j) ~= cut(dimx,j+1)
        line([j+0.5+padding j+0.5+padding],[dimx-0.5+padding dimx+0.5+padding],'LineWidth',myWidth,'Color',squeeze(colorList(1,:)));
    end
end
for i=1:dimx %fill in edges
    if cut(i,1)
        line([0.5 0.5]+padding,[i-0.5 i+0.5]+padding,'LineWidth',myWidth);
    end
    if cut(i,dimy)
        line([dimy+0.5 dimy+0.5]+padding,[i-0.5 i+0.5]+padding,'LineWidth',myWidth);
    end                            
end
for j=1:dimy
    if cut(1,j)
        line([j-0.5 j+0.5]+padding,[0.5 0.5]+padding,'LineWidth',myWidth);
    end
    if cut(dimx,j)
        line([j-0.5 j+0.5]+padding,[dimx+0.5 dimx+0.5]+padding,'LineWidth',myWidth);
    end
end
hold on; plot(queryPoint(2)+padding,queryPoint(1)+padding,'bp','MarkerSize',15,'MarkerFaceColor','b');
%print('-painters','-dpng','-r300','~/matlab_output/figure_specified_point.png');
%print('-painters','-dpdf','-r300','~/matlab_output/figure_specified_point.pdf');

figure; semilogx(emdList(:,3),emdList(:,4),'b-x','MarkerSize',10,'LineWidth',1);
xlim([0.5 50]);
xlabel('Region scale (mm)');
ylabel('Signal intensity heterogeneity (EMD)');
title('Heterogeneity vs. scale at a point');
figureHandle=gcf;
set(findall(figureHandle,'type','text'),'fontSize',16);
set(gca,'FontSize',16);
%print('-painters','-dpdf','-r300','~/matlab_output/figure_specified_point_graph.pdf');

