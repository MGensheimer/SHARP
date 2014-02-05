function [partsCell,partsCellLen,parentCell]=partition2D(myImage,myMask,sigmai,sigmax_constant,voxelSize)

% 2D SHARP partitioning algorithm
%
% Inputs:
% myImage: 2D double-valued image
% myMask: 2D image. Values of 1 inside tumor, 0 outside tumor.
% sigmai: intensity sensitivity constant for Normalized Cut
% sigmax_constant: distance sensitivity constant for Normalized Cut
% voxelSize: 1x2 array of voxel dimensions in mm

% Outputs:
% partsCell: If there are n levels of divisions, this is a
% [n+1]-dimensional cell array. Each cell holds a 2-dimensional array. Each
% row in this array lists the voxels included in a partition. Smaller
% partitions are padded with 0's.
% partsCellLen: [n+1]-dimensional cell array. Each cell holds a vector with
% the number of voxels in each partition at that division level.
% parentCell: [n+1]-dimensional cell array. This is a "poor man's
% tree" data structure which can be used to reconstruct the SHARP partition
% tree. It has parallel structure to partsCell. To find the parent
% partition of a partition in partsCell, go to the corresponding location
% in parentCell and look at the number index. This number indicates the row
% index of the parent partition in partsCell (at 1 lower level).

rng(0); %give same seed to kmeans and randperm each time algorithm is run
[dimx, dimy]=size(myImage);
whereTumor=find(myMask>0);
whereTumorInv=zeros(dimx*dimy,1,'uint32');
whereTumorInv(whereTumor)=1:size(whereTumor,1);

[whereTumorX,whereTumorY]=find(myMask>0);
myFilter=[0 1 0; 1 1 1; 0 1 0];

warning('off','stats:kmeans:EmptyCluster');
warning('off','stats:kmeans:EmptyClusterInBatchUpdate');
numSamples=300;
numKmeansReps=10;
voxelSize=double(voxelSize);
needsDivide=zeros(ceil(size(whereTumor,1)/3),1,'uint32'); needsDivide(1)=1;
needsDivideCount=1;
level=1;
maxLevel=50;
partsCell=cell(maxLevel,1);
partsCellLen=cell(maxLevel,1);
parentCell=cell(maxLevel,1);
partsCell{1}=1:size(whereTumor,1);
partsCellLen{1}=size(whereTumor,1);
parentCell{1}=0;
while needsDivideCount && level <= maxLevel-1
    level=level+1;
    %max # partitions in this level is 1/2 of prior level; max partition
    %size in this level is biggest size in prior level - 1
    levelParts=zeros(size(partsCell{level-1},1)*2,max(partsCellLen{level-1})-1,'uint32');
    levelPartsLen=zeros(size(partsCell{level-1},1)*2,1,'uint32');
    toDivide=needsDivide;
    toDivideCount=needsDivideCount;
    needsDivide=zeros(ceil(size(whereTumor,1)/3),1,'uint32');
    needsDivideCount=0;
    for currentPart=1:toDivideCount
        labelsIn=partsCell{level-1}(toDivide(currentPart),1:partsCellLen{level-1}(toDivide(currentPart)));
        sigmax=sigmax_constant*sqrt(size(labelsIn,2)*prod(voxelSize));
        perm=randperm(size(labelsIn,2));
        [~,reversePerm]=sort(perm);
        if size(labelsIn,2) == 2
            for child=1:2
                levelParts((currentPart-1)*2+child,1)=labelsIn(child);
                levelPartsLen((currentPart-1)*2+child)=1;
                parentCell{level}((currentPart-1)*2+child)=toDivide(currentPart);
            end
        else
            if size(labelsIn,2) <= numSamples
                [x1, x2]=meshgrid(whereTumorX(labelsIn(perm)));
                [y1, y2]=meshgrid(whereTumorY(labelsIn(perm)));
                [int1, int2]=meshgrid(myImage(whereTumor(labelsIn(perm))));
                distMat=exp(-( ( ((x2-x1)*voxelSize(1)).^2 + ((y2-y1)*voxelSize(2)).^2 )/(sigmax^2) +((int1-int2)/sigmai).^2));
                [V,~]=ncut(distMat);
            else
                [x1, x2]=meshgrid(whereTumorX(labelsIn(perm)), whereTumorX(labelsIn(perm(1:numSamples))));
                [y1, y2]=meshgrid(whereTumorY(labelsIn(perm)), whereTumorY(labelsIn(perm(1:numSamples))));
                [int1, int2]=meshgrid(myImage(whereTumor(labelsIn(perm))), myImage(whereTumor(labelsIn(perm(1:numSamples)))));
                distMat=exp(-( ( ((x2-x1)*voxelSize(1)).^2 + ((y2-y1)*voxelSize(2)).^2 )/(sigmax^2) +((int1-int2)/sigmai).^2));
                [V,~]=nystrom_ncut(distMat(1:numSamples,1:numSamples),distMat(1:numSamples,numSamples+1:end));
            end
            V=V(reversePerm,:);
            V1=V(:,1:2); % choose the 2 eigenvectors
            Vnormalized=V1./repmat(sqrt(sum(V1.^2,2)),1,2);
            Vnormalized(isnan(Vnormalized))=0;
            Vnormalized(isinf(Vnormalized))=0;
            labels=kmeans(Vnormalized,2,'replicates',numKmeansReps);
            if size(unique(labels),1)==1 % check if one of the 2 partitions is empty (bad!)
                keyboard;
            end
            myCut=zeros(dimx,dimy); % find # of connected components (want to be 2)
            myCut(whereTumor(labelsIn))=labels;
            myCutExpanded=zeros(dimx+8,dimy+8);
            myCutExpanded(5:end-4,5:end-4)=myCut;
            CC1=bwconncomp(myCutExpanded==1,4);
            CC2=bwconncomp(myCutExpanded==2,4);
            numObjects = [CC1.NumObjects CC2.NumObjects];
            if CC1.NumObjects+CC2.NumObjects > 2 % one partition is discontinuous
                for whichPart=1:2
                    if numObjects(whichPart) > 1
                        if whichPart==1
                            CC=CC1;
                        else
                            CC=CC2;
                        end
                        labels=labelmatrix(CC);
                        otherMask=myCutExpanded==3-whichPart;
                        numPixels = cellfun(@numel,CC.PixelIdxList);
                        [~,largestPartID]=max(numPixels);
                        for whichObject=1:CC.NumObjects
                            if whichObject ~= largestPartID
                                objectExpanded=imfilter(labels==whichObject,myFilter);
                                if max(max(objectExpanded & otherMask))
                                    myCutExpanded(CC.PixelIdxList{whichObject})=3-whichPart; %assign object to other partition
                                end
                            end
                        end
                    end
                end
                bestCut=myCutExpanded(5:end-4,5:end-4);
                part{1}=whereTumorInv(bestCut==1)';
                part{2}=whereTumorInv(bestCut==2)';
            else
                part{1}=labelsIn(labels == 1);
                part{2}=labelsIn(labels == 2);                    
            end
            for child=1:2
                levelParts((currentPart-1)*2+child,1:size(part{child},2))=part{child};
                levelPartsLen((currentPart-1)*2+child)=size(part{child},2);
                parentCell{level}((currentPart-1)*2+child)=toDivide(currentPart);
                if size(part{child},2)>1
                    needsDivideCount=needsDivideCount+1;
                    needsDivide(needsDivideCount)=(currentPart-1)*2+child;
                end
            end
        end
    end
    [u,v,~]=find(levelParts);
    partsCell{level}=levelParts(1:max(u),1:max(v));
    partsCellLen{level}=levelPartsLen(1:max(u));
end
[temp,~]=find(cellfun(@numel,partsCell));
partsCell=partsCell(1:max(temp));
partsCellLen=partsCellLen(1:max(temp));
parentCell=parentCell(1:max(temp));