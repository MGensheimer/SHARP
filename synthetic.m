%Validates SHARP algorithm using synthetic data: three classes of images
%with varying scales of heterogeneity

voxelSize=[1 1];
dimx=64; dimy=64;
rect_levels=12;

% output example images for the 3 classes
figure;
intensities=2.^(rect_levels-1:-1:0);
myImage=zeros(dimx,dimy);
for level=1:rect_levels
    x_rect_size=dimx/(2^ceil(level/2));
    y_rect_size=dimy/(2^ceil((level-1)/2));
    for i=0:(2^ceil(level/2))-1
        for j=0:(2^ceil((level-1)/2))-1
            myImage(1+i*x_rect_size:(i+1)*x_rect_size,1+j*y_rect_size:(j+1)*y_rect_size) = ...
            myImage(1+i*x_rect_size:(i+1)*x_rect_size,1+j*y_rect_size:(j+1)*y_rect_size) + (rand-0.5)*intensities(level);
        end
    end
end
subplot(3,1,1); imshow(myImage,[]); title(gca,'Large blocks');

intensities=ones(1,rect_levels);
myImage=zeros(dimx,dimy);
for level=1:rect_levels
    x_rect_size=dimx/(2^ceil(level/2));
    y_rect_size=dimy/(2^ceil((level-1)/2));
    for i=0:(2^ceil(level/2))-1
        for j=0:(2^ceil((level-1)/2))-1
            %[level x_rect_size y_rect_size i j]
            myImage(1+i*x_rect_size:(i+1)*x_rect_size,1+j*y_rect_size:(j+1)*y_rect_size) = ...
            myImage(1+i*x_rect_size:(i+1)*x_rect_size,1+j*y_rect_size:(j+1)*y_rect_size) + (rand-0.5)*intensities(level);
        end
    end
end
subplot(3,1,2); imshow(myImage,[]); title(gca,'Large+small blocks');

intensities=2.^(0:rect_levels-1);
myImage=zeros(dimx,dimy);
for level=1:rect_levels
    x_rect_size=dimx/(2^ceil(level/2));
    y_rect_size=dimy/(2^ceil((level-1)/2));
    for i=0:(2^ceil(level/2))-1
        for j=0:(2^ceil((level-1)/2))-1
            %[level x_rect_size y_rect_size i j]
            myImage(1+i*x_rect_size:(i+1)*x_rect_size,1+j*y_rect_size:(j+1)*y_rect_size) = ...
            myImage(1+i*x_rect_size:(i+1)*x_rect_size,1+j*y_rect_size:(j+1)*y_rect_size) + (rand-0.5)*intensities(level);
        end
    end
end
subplot(3,1,3); imshow(myImage,[],'Border','tight'); title(gca,'Small blocks');

%create plot of EMD versus scale for 5 example images of each class
patterns=3;
reps=5;
intensities=zeros(patterns,rect_levels);
intensities(1,:)=2.^(rect_levels-1:-1:0);
intensities(2,:)=ones(1,rect_levels);
intensities(3,:)=2.^(0:rect_levels-1);
sigmai_constant=2.5;
sigmax_constant=0.5;
numBinnedLevels=10;
meanyAllClasses=zeros(patterns,numBinnedLevels);
stdyAllClasses=zeros(patterns,numBinnedLevels);
for whichIntensities=1:patterns
    yall=zeros(reps,numBinnedLevels);
    for rep=1:reps
        rep
        myImage=zeros(dimx,dimy);
        for level=1:rect_levels
            x_rect_size=dimx/(2^ceil(level/2));
            y_rect_size=dimy/(2^ceil((level-1)/2));
            for i=0:(2^ceil(level/2))-1
                for j=0:(2^ceil((level-1)/2))-1
                    myImage(1+i*x_rect_size:(i+1)*x_rect_size,1+j*y_rect_size:(j+1)*y_rect_size) = ...
                    myImage(1+i*x_rect_size:(i+1)*x_rect_size,1+j*y_rect_size:(j+1)*y_rect_size) + (rand-0.5)*intensities(whichIntensities,level);
                end
            end
        end
        myMask=ones(dimx,dimy);
        whereTumor=find(myMask>0);
        tumorIntensity=myImage(whereTumor);
        sortedIntensity=sort(tumorIntensity);
        quart1Intensity=sortedIntensity(ceil(size(sortedIntensity,1)*0.25));
        quart3Intensity=sortedIntensity(ceil(size(sortedIntensity,1)*0.75));
        sigmai=(quart3Intensity-quart1Intensity)*sigmai_constant/1.4; %interquartile / 1.4 = std dev for normal dist
        [partsCell,partsCellLen,parentCell]=partition2D(myImage,myMask,sigmai,sigmax_constant,voxelSize);
        [whereTumorX,whereTumorY]=find(myMask>0);
        numLevels=size(partsCell,1)-1;
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
                yall(rep,i)=y(i);
            end
        end
    end
    meanyall=mean(yall,1);
    stdyall=std(yall,1);
    meanyAllClasses(whichIntensities,:)=meanyall;
    stdyAllClasses(whichIntensities,:)=stdyall;
end
figure;
h=errorbar(log2(repmat(x(x>0)',patterns,1)'),meanyAllClasses(:,x>0)',stdyAllClasses(:,x>0)','-o','LineWidth',1);%,'-o','MarkerFaceColor','r');
set(gca,'XTickLabel',[]);
xt=get(gca,'XTick');
y1=get(gca,'YLim');
str=cellstr(num2str(xt(:),'2^{%d}'));
hTxt=text(xt,y1(ones(size(xt))),str,'Interpreter','tex','VerticalAlignment','top','HorizontalAlignment','center');
xlabel('Region scale (mm)');
xlabh=get(gca,'XLabel'); set(xlabh,'Position',get(xlabh,'Position')-[0 .3 0]);
ylabel('Signal intensity heterogeneity (mean EMD)');
legend('Large blocks','Large+small blocks','Small blocks','location','NorthWest');
