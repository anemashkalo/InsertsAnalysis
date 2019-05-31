function [mask2,reader] = get_mask(dir_img,chan1,pos,tgroup,thresh,foreground,simple_seg)

ff = readAndorDirectory(dir_img);
tpname = getAndorFileName(ff,pos,tgroup,[],ff.w(chan1));%chan
if simple_seg
filename = [tpname(1:end-4), '_Simple Segmentation.h5'];%
else
filename = [tpname(1:end-4), '_Probabilities.h5'];%Probabilities
end
if exist(filename)
reader = bfGetReader(tpname);
nz=reader.getSizeZ;
nT = reader.getSizeT;
nz=reader.getSizeZ;
% for Simple Segmentation input
if ~contains(filename,'Probabilities')
immask = h5read(filename, '/exported_data');
immask2 = squeeze(immask(1,:,:,:));
mask2 = zeros(size(immask2,1),size(immask2,2),size(immask2,3));
clear mask2
close all
for jj=1:size(immask2,3)
t = ((immask2(:,:,jj)'));
t(t==1)=0;
t(t==2)=1;
%imshow(mask2(:,:,jj),[]);
mask2(:,:,jj) = im2bw(t,0);
end
else
% for Probability masks
[mask2]=readIlastikProbMask(filename,thresh,foreground);
end
else %only return the image
    disp('Ilastik files were not found, returning only image reader')
    mask2 =[];
  reader = bfGetReader(tpname);
nz=reader.getSizeZ;
nT = reader.getSizeT;
nz=reader.getSizeZ;

    
end



end


