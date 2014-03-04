info=dicominfo('mri_image.dcm');
voxelSize=zeros(1,2);
voxelSize(1:2)=info.PixelSpacing;
imageUncropped=double(dicomread('mri_image.dcm'));
maskUncropped=double(dicomread('mri_mask.dcm'));
[row,col]=find(maskUncropped>0);
image=double(imageUncropped(min(row):max(row),min(col):max(col)));
mask=double(maskUncropped(min(row):max(row),min(col):max(col)));
mask(mask<0)=0;
[dimx dimy]=size(image);
image(~mask)=NaN;
warning('off','images:graycomatrix:scaledImageContainsNan');

histLevels=32; %used in Gibbs 2003
% calculate GLCM for angles 0, 45, 90, 135
glcm1=graycomatrix(image,'NumLevels',histLevels,'GrayLimits',[],'Offset',[0 1],'Symmetric',true);
glcm2=graycomatrix(image,'NumLevels',histLevels,'GrayLimits',[],'Offset',[-1 1],'Symmetric',true);
glcm3=graycomatrix(image,'NumLevels',histLevels,'GrayLimits',[],'Offset',[-1 0],'Symmetric',true);
glcm4=graycomatrix(image,'NumLevels',histLevels,'GrayLimits',[],'Offset',[-1 -1],'Symmetric',true);
glcmAll=cat(3,glcm1,glcm2,glcm3,glcm4);
glcmMean=mean(glcmAll,3); %find average of the four GLCMs
stats=GLCM_Features4(glcmMean,0); %pairs flag set to 0 because Symmetric is true
stats.entro %entropy
stats.senth %sum entropy
stats.sosvh %variance
