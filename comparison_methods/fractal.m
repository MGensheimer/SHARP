fractal_levels=4;

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
image(~mask)=0;

dim_start=ceil(max(dimx,dimy)/2^(fractal_levels-1))*2^(fractal_levels-1);
image2=zeros(dim_start,dim_start);
image2(1:dimx,1:dimy)=image;
voxres=[]; cvres=[];
for level=1:fractal_levels
    newImage=imresize(image2,0.5^(level-1),'bilinear','Antialiasing',false);
    voxres(level)=size(newImage(:),1);
    cvres(level)=std(newImage(:))/mean(newImage(:));
end
slope=log(cvres(1)./cvres)' \ log(voxres/voxres(1))';
fd=1-(1/slope);
ff=sum((log(voxres/voxres(1))-log(cvres(1)./cvres)*slope).^2);
fd %fractal dimension
ff %fractal fit